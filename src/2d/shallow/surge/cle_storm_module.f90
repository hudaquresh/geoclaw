! ==============================================================================
! storm_module 
!
! Module contains routines for constructing a wind and pressure field based on
! the CLE hurricane module. 
!
! Many routines are based loosely on PADCIRC version 45.12 03/17/2006, but 
! updated for the current CLE model. 
! ==============================================================================
module cle_storm_module

    implicit none 
    save
    
    ! Chavas storm type definition    
    type cle_storm_type 
        ! Fore/hindcast size and current position 
        integer :: num_casts 

        ! These paramters are located at time points but are interpolated in  
        ! time and space when the relevant fields are requested. 

        ! Location of storm 
        ! Track is a triplet with (time,longitude,latitude)  
        real(kind=8), allocatable :: track(:,:) 

        ! Storm physics 
        real(kind=8), allocatable :: max_wind_radius(:)
        real(kind=8), allocatable :: max_wind_speed(:)
        real(kind=8), allocatable :: central_pressure(:)
        real(kind=8), allocatable :: ambient_pressure(:)
        real(kind=8), allocatable :: rrp(:) 

        ! Approximate velocity of storm, approximated via the track points
        ! using a first order difference on the sphere 
        real(kind=8), allocatable :: velocity(:,:) 

    end type cle_storm_type

    logical, private :: module_setup = .false.  

    ! Internal tracking variables for storm 
    integer, private :: last_storm_index  

    logical, private :: DEBUG = .false.

    ! Atmospheric boundary layer, input variable in ADCIRC but always is 
    ! set to the following value 
    real(kind=8), parameter :: atmos_boundary_layer = 0.9d0 

    ! Sampling adjustment from 1 min to 10 min winds 
    real(kind=8), parameter :: sampling_time = 0.88d0 

    ! Storm field ramping width - Represents curdely the ramping radial area 
    real(kind=8), parameter :: RAMP_WIDTH = 100.0d3  

    ! Global constants #TODO: Some of these are in geoclaw already  
    real(kind=8) :: pi=3.1415927 
    real(kind=8) :: omega=7.2921e-5
    real(kind=8) :: chi=1.0
    real(kind=8) :: alpha=1.0
 
contains

    ! Setup routine for the cle model 
    subroutine set_cle_storm(storm_data_path, storm, log_unit)

        use geoclaw_module, only: deg2rad, spherical_distance, coordinate_system
        use amr_module, only: t0, rinfinity

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(cle_storm_type), intent(in out) :: storm
        integer, intent(in) :: log_unit

        ! Local storage
        integer, parameter :: data_file = 701
        integer :: i, k, io_status, num_casts
        real(kind=8) :: forecast_time,last_time,x(2),y(2),ds,dt,dx,dy,theta

        ! Reading buffer variables
        integer :: year,month,day,hour,forecast,lat,lon,max_wind_speed
        integer :: central_pressure,RRP,max_wind_radius
        character(len=4) :: cast_type
        character(len=1) :: direction(2)

        ! Note that the JAM file format has not been tested yet and will later
        ! be added as an option.
        character(len=4), parameter :: file_format = "NOAA"

        ! File format string
        character(len=*), parameter :: JMA_FORMAT = "(i2,i2,i2,i2,8x,i3,1x,"//&
                            "i4,1x,i4,5x,i3)"
        character(len=*), parameter :: NOAA_FORMAT = "(8x,i4,i2,i2,i2,6x,a4,"//&
                            "2x,i3,1x,i4,a1,2x,i4,a1,2x,i3,2x,i4,47x,i3,2x,i3)"

        if (.not. module_setup) then
            
            ! Storm type only works on lat-long coordinate systems
            if (coordinate_system /= 2) then
                stop "Holland storm type does only works on lat-long coordinates."
            endif

            ! Open data file
            if (present(storm_data_path)) then
                print *,'Reading storm date file ',storm_data_path
                open(unit=data_file,file=storm_data_path,status='old', &
                     action='read',iostat=io_status)
            else
                print *,'Reading storm date file ./storm.data'
                open(unit=data_file,file="./storm.data",status='old', &
                     action='read',iostat=io_status)
            endif
            if (io_status /= 0) then
                print "(a,i2)", "Error opening storm data file. status = ", io_status
                stop 
            endif            

            ! Count number of data lines
            num_casts = 0
            last_time = -rinfinity
            do
                if (file_format == "NOAA") then
                    read(data_file,fmt=NOAA_FORMAT,iostat=io_status) year,month,day, &
                        hour,cast_type,forecast,lat,direction(2),lon,direction(1), &
                        max_wind_speed,central_pressure,RRP,max_wind_radius
                else if (file_format == "JAM") then
                    ! JAM may be missing RRP parameter, may need to set this based
                    ! on other data in the file.  It is only used in the field 
                    ! ramping function so it might not be an issue
                    read(data_file,fmt=JMA_FORMAT,iostat=io_status) year, month, day, &
                            hour, lat, lon, central_pressure, max_wind_speed, max_wind_radius
                else
                    print *,"ERROR - Unrecognized storm data file format."
                    stop
                endif

                ! Exit loop if we ran into an error or we reached the end of the file
                if (io_status /= 0) exit

                ! Skip counting this line if time is repeated
                forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
                if (abs(forecast_time - last_time) >= 1.8d3) then
                    num_casts = num_casts + 1
                endif
                last_time = forecast_time
            end do
            rewind(data_file)

            write(log_unit,"('Forecasts = ',i3)") num_casts

            ! Allocate storm parameter file variables
            allocate(storm%track(3,num_casts))
            allocate(storm%max_wind_speed(num_casts))
            allocate(storm%max_wind_radius(num_casts))
            allocate(storm%central_pressure(num_casts))
            allocate(storm%rrp(num_casts))

            ! Now re-read the file's contents
            i = 0
            do while (i < num_casts)
                if (file_format == "NOAA") then
                    read(data_file,fmt=NOAA_FORMAT) year,month,day,hour,cast_type, &
                        forecast,lat,direction(2),lon,direction(1),max_wind_speed, &
                        central_pressure,RRP,max_wind_radius
                else if (file_format == "JAM") then
                    read(data_file,fmt=JMA_FORMAT,iostat=io_status) year, month, day, &
                            hour, lat, lon, central_pressure, max_wind_speed, max_wind_radius
                else
                    print *,"ERROR - Unrecognized storm data file format."
                    stop
                endif


                ! Skip counting this line if time is repeated
                forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
                if (abs(forecast_time - last_time) < 1.8d3) then
                    cycle
                endif
                i = i + 1
                last_time = forecast_time

                ! Storm position
                ! Conversions:
                !  lon - Convert 10ths of degs to degs, depends on E,W
                !  lat - Convert 10ths of degs to degs
                storm%track(1,i) = date_to_seconds(year,month,day,hour,0,0.d0)
                if (direction(1) == "E") then
                    storm%track(2,i) = real(lon,kind=8) / 10.d0
                else
                    storm%track(2,i) = -real(lon,kind=8) / 10.d0
                endif
                if (direction(2) == "N") then
                    storm%track(3,i) = real(lat,kind=8) / 10.d0
                else
                    storm%track(3,i) = -real(lat,kind=8) / 10.d0
                endif

                ! Storm intensity
                ! Conversions:
                !  max_wind_speed - Convert knots to m/s
                !  max_wind_radius  - convert from nm to m
                !  central_pressure - convert from mbar to Pa
                !  Radius of last isobar contour - convert from nm to m
                storm%max_wind_speed(i) = real(max_wind_speed,kind=8) * 0.51444444d0
                storm%max_wind_radius(i) = real(max_wind_radius,kind=8) * 1.852000003180799d0 * 1000.d0
                storm%central_pressure(i) = real(central_pressure,kind=8) * 100.d0
                storm%rrp(i) = real(RRP,kind=8) * 1.852000003180799d0 * 1000.d0

            enddo

            ! Calculate storm speed 
            allocate(storm%velocity(2,num_casts))
            do i=1,num_casts - 1
                ! Calculate velocity based on great circle distance between

                ! locations of storm
                x = storm%track(2:3,i)
                y = storm%track(2:3,i+1)

                dt = storm%track(1,i + 1) - storm%track(1,i)

                ds = spherical_distance(x(1), 0.5d0 * (x(2) + y(2)), &
                                        y(1), 0.5d0 * (x(2) + y(2)))
                storm%velocity(1,i) = sign(ds / dt,y(1) - x(1))

                
                ds = spherical_distance(0.5d0 * (x(1) + y(1)), x(2), &
                                        0.5d0 * (x(1) + y(1)), y(2))
                storm%velocity(2,i) = sign(ds / dt,y(2) - x(2))
            end do

            ! Use last approximation for velocity point going forward
            storm%velocity(:,num_casts) = storm%velocity(:,num_casts - 1)

            ! Record number of casts
            storm%num_casts = num_casts

            if (t0 < storm%track(1,1)) then
                print *,t0,storm%track(1,1)
                stop "Start time is before first forecast time."
            endif

            ! This is used to speed up searching for correct storm data
            last_storm_index = 2
            last_storm_index = storm_index(t0,storm)
            if (last_storm_index == -1) then
                print *,"Forecast not found for time ",t0,'.'
                stop
            endif

            ! Log everything to the surge log file
            write(log_unit,*) ""
            write(log_unit,*) "Storm Track and Strength"
            write(log_unit,*) ""
            do i=1,storm%num_casts
                write(log_unit,"(8e26.16)") (storm%track(k,i),k=1,3),  &
                                            (storm%velocity(k,i),k=1,2), &
                                             storm%max_wind_radius(i), &
                                             storm%max_wind_speed(i),  &
                                             storm%central_pressure(i)
            enddo

            module_setup = .true.
        end if
!
!        ! Storm type only works on lat-long coordinate systems
!        if (coordinate_system /= 2) then
!            stop "CLE storm type does only works on lat-long coordinates."
!        endif
!
!        ! Open data file
!        if (present(storm_data_path)) then
!            print *,'Reading storm date file ',storm_data_path
!            open(unit=data_file,file=storm_data_path,status='old', &
!                 action='read',iostat=io_status)
!        else
!            print *,'Reading storm date file ./storm.data'
!            open(unit=data_file,file="./storm.data",status='old', &
!                 action='read',iostat=io_status)
!        endif
!        if (io_status /= 0) then
!            print "(a,i2)", "Error opening storm data file. status = ", io_status
!            stop 
!        endif            
!
!        ! Count number of data lines
!        num_casts = 0
!        last_time = -rinfinity
!        do
!            if (file_format == "NOAA") then
!                read(data_file,fmt=NOAA_FORMAT,iostat=io_status) year,month,day, &
!                    hour,cast_type,forecast,lat,direction(2),lon,direction(1), &
!                    max_wind_speed,central_pressure,RRP,max_wind_radius
!            else if (file_format == "JAM") then
!                ! JAM may be missing RRP parameter, may need to set this based
!                ! on other data in the file.  It is only used in the field 
!                ! ramping function so it might not be an issue
!                read(data_file,fmt=JMA_FORMAT,iostat=io_status) year, month, day, &
!                        hour, lat, lon, central_pressure, max_wind_speed, max_wind_radius
!            else
!                print *,"ERROR - Unrecognized storm data file format."
!                stop
!            endif
!
!            ! Exit loop if we ran into an error or we reached the end of the file
!            if (io_status /= 0) exit
!
!            ! Skip counting this line if time is repeated
!            forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
!            if (abs(forecast_time - last_time) >= 1.8d3) then
!                num_casts = num_casts + 1
!            endif
!            last_time = forecast_time
!        end do
!        rewind(data_file)
!
!        write(log_unit,"('Forecasts = ',i3)") num_casts
!
!        ! Allocate storm parameter file variables
!        allocate(storm%track(3,num_casts))
!        allocate(storm%max_wind_speed(num_casts))
!        allocate(storm%max_wind_radius(num_casts))
!        allocate(storm%central_pressure(num_casts))
!        allocate(storm%rrp(num_casts))
!
!        ! Now re-read the file's contents
!        i = 0
!        do while (i < num_casts)
!            if (file_format == "NOAA") then
!                read(data_file,fmt=NOAA_FORMAT) year,month,day,hour,cast_type, &
!                    forecast,lat,direction(2),lon,direction(1),max_wind_speed, &
!                    central_pressure,RRP,max_wind_radius
!            else if (file_format == "JAM") then
!                read(data_file,fmt=JMA_FORMAT,iostat=io_status) year, month, day, &
!                        hour, lat, lon, central_pressure, max_wind_speed, max_wind_radius
!            else
!                print *,"ERROR - Unrecognized storm data file format."
!                stop
!            endif
!
!
!            ! Skip counting this line if time is repeated
!            forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
!            if (abs(forecast_time - last_time) < 1.8d3) then
!                cycle
!            endif
!            i = i + 1
!            last_time = forecast_time
!
!            ! Storm position
!            ! Conversions:
!            !  lon - Convert 10ths of degs to degs, depends on E,W
!            !  lat - Convert 10ths of degs to degs
!            storm%track(1,i) = date_to_seconds(year,month,day,hour,0,0.d0)
!            if (direction(1) == "E") then
!                storm%track(2,i) = real(lon,kind=8) / 10.d0
!            else
!                storm%track(2,i) = -real(lon,kind=8) / 10.d0
!            endif
!            if (direction(2) == "N") then
!                storm%track(3,i) = real(lat,kind=8) / 10.d0
!            else
!                storm%track(3,i) = -real(lat,kind=8) / 10.d0
!            endif
!
!            ! Storm intensity
!            ! Conversions:
!            !  max_wind_speed - Convert knots to m/s
!            !  max_wind_radius  - convert from nm to m
!            !  central_pressure - convert from mbar to Pa
!            !  Radius of last isobar contour - convert from nm to m
!            storm%max_wind_speed(i) = real(max_wind_speed,kind=8) * 0.51444444d0
!            storm%max_wind_radius(i) = real(max_wind_radius,kind=8) * 1.852000003180799d0 * 1000.d0
!            storm%central_pressure(i) = real(central_pressure,kind=8) * 100.d0
!            storm%rrp(i) = real(RRP,kind=8) * 1.852000003180799d0 * 1000.d0
!
!        enddo
!
!        ! Calculate storm speed 
!        allocate(storm%velocity(2,num_casts))
!        do i=1,num_casts - 1
!            ! Calculate velocity based on great circle distance between
!
!            ! locations of storm
!            x = storm%track(2:3,i)
!            y = storm%track(2:3,i+1)
!
!            dt = storm%track(1,i + 1) - storm%track(1,i)
!
!            ds = spherical_distance(x(1), 0.5d0 * (x(2) + y(2)), &
!                                    y(1), 0.5d0 * (x(2) + y(2)))
!            storm%velocity(1,i) = sign(ds / dt,y(1) - x(1))
!
!            
!            ds = spherical_distance(0.5d0 * (x(1) + y(1)), x(2), &
!                                    0.5d0 * (x(1) + y(1)), y(2))
!            storm%velocity(2,i) = sign(ds / dt,y(2) - x(2))
!        end do
!
!        ! Use last approximation for velocity point going forward
!        storm%velocity(:,num_casts) = storm%velocity(:,num_casts - 1)
!
!        ! Record number of casts
!        storm%num_casts = num_casts
!
!        if (t0 < storm%track(1,1)) then
!            print *,t0,storm%track(1,1)
!            stop "Start time is before first forecast time."
!        endif
!
!        ! This is used to speed up searching for correct storm data
!        last_storm_index = 2
!        last_storm_index = storm_index(t0,storm)
!        if (last_storm_index == -1) then
!            print *,"Forecast not found for time ",t0,'.'
!            stop
!        endif
!
!        ! Log everything to the surge log file
!        write(log_unit,*) ""
!        write(log_unit,*) "Storm Track and Strength"
!        write(log_unit,*) ""
!        do i=1,storm%num_casts
!            write(log_unit,"(8e26.16)") (storm%track(k,i),k=1,3),  &
!                                        (storm%velocity(k,i),k=1,2), &
!                                         storm%max_wind_radius(i), &
!                                         storm%max_wind_speed(i),  &
!                                         storm%central_pressure(i)
!        enddo
!       
!        module_setup = .true. 
!     end if 

    end subroutine set_cle_storm


    ! ==========================================================================
    !  real(kind=8) pure date_to_seconds(year,months,days,hours,minutes,seconds)
    !    Convert time from year, month, day, hour, min, sec to seconds since the
    !    beginning of the year.
    ! ==========================================================================
    pure real(kind=8) function date_to_seconds(year,months,days,hours,minutes, &
                                               seconds) result(time)
      
        implicit none

        ! Input
        integer, intent(in) :: year, months, days, hours, minutes
        real(kind=8), intent(in) :: seconds

        ! Local storage
        integer :: total_days

        ! Count number of days
        total_days = days

        ! Add days for months that have already passed
        if (months > 1) total_days = total_days + 31
        if (months > 2) then
            if (int(year / 4) * 4 == year) then
                total_days = total_days + 29
            else
                total_days = total_days + 28
            endif
        endif
        if (months > 3)  total_days = total_days + 31
        if (months > 4)  total_days = total_days + 30
        if (months > 5)  total_days = total_days + 31
        if (months > 6)  total_days = total_days + 30
        if (months > 7)  total_days = total_days + 31
        if (months > 8)  total_days = total_days + 31
        if (months > 9)  total_days = total_days + 30
        if (months > 10) total_days = total_days + 31
        if (months > 11) total_days = total_days + 30

        ! Convert everything to seconds since the beginning of the year
        time = real((total_days - 1) * 86400 + hours * 3600 + minutes * 60,kind=8)
        time = time + seconds

    end function date_to_seconds


    ! ==========================================================================
    !  cle_storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function cle_storm_location(t,storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(cle_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        ! Junk storage
        real(kind=8) :: junk(2)

        call get_cle_storm_data(t,storm,location, &
                                        junk,junk(1),junk(1),junk(1),junk(1))

    end function cle_storm_location

    ! ==========================================================================
    !  cle_storm_direction
    !   Angle off of due north that the storm is traveling
    ! ==========================================================================
    real(kind=8) function cle_storm_direction(t, storm) result(theta)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(cle_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: junk(2), velocity(2)

        ! Fetch velocity of storm which has direction encoded in it
        call get_cle_storm_data(t, storm, junk, velocity, junk(1),  &
                                                    junk(1), junk(1), junk(1))

        ! Unit directional vector
        theta = atan2(velocity(2),velocity(1))

    end function cle_storm_direction

    ! ==========================================================================
    !  storm_index(t,storm)
    !    Finds the index of the next storm data point
    ! ==========================================================================
    integer pure function storm_index(t,storm) result(index)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(cle_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: t0,t1
        logical :: found

        ! Figure out where we are relative to the last time we checked for the
        ! index (stored in last_storm_index)
        
        ! Check if we are already beyond the end of the last forecast time
        if (last_storm_index == storm%num_casts + 1) then
            index = storm%num_casts + 1
        else
            t0 = storm%track(1,last_storm_index - 1)
            t1 = storm%track(1,last_storm_index)
            if (t0 < t .and. t <= t1) then
                index = last_storm_index
            else if ( t1 < t ) then
                found = .false.
                do index=last_storm_index+1,storm%num_casts
                    if (t < storm%track(1,index)) then
                        found = .true.
                        exit
                    endif
                enddo
                ! Assume we have gone past last forecast time
                if (.not. found) then
                    index = storm%num_casts + 1
                endif
            else ! t <= t0
                if (last_storm_index == 2) then
                    index = -1
                else
                    do index=last_storm_index-1,2,-1
                        if (storm%track(1,index-1) < t) exit
                    enddo
                endif
            endif
        endif

    end function storm_index
    
    ! ==========================================================================
    !  get_cle_storm_data()
    !    Interpolates in time and returns storm data.
    ! ==========================================================================
    subroutine get_cle_storm_data(t, storm, location, velocity,   &
                                                max_wind_radius,  &
                                                max_wind_speed,   &
                                                central_pressure, &
                                                rrp)

        use geoclaw_module, only: deg2rad, latlon2xy, xy2latlon

        implicit none

        ! Input
        real(kind=8), intent(in) :: t                       ! Current time
        type(cle_storm_type), intent(in) :: storm   ! Storm

        ! Output
        real(kind=8), intent(out) :: location(2), velocity(2)
        real(kind=8), intent(out) :: max_wind_radius, max_wind_speed
        real(kind=8), intent(out) :: central_pressure, rrp

        ! Local
        real(kind=8) :: fn(8), fnm(8), weight, tn, tnm, x(2)
        integer :: i

        ! Increment storm data index if needed and not at end of forecast
        i = storm_index(t,storm)
        last_storm_index = i

        ! List of possible error conditions
        if (i <= 1) then
            if (i == 0) then        
                print *,"Invalid storm forecast requested for t = ",t
                print *,"Time requested is before any forecast data."
                print *,"    first time = ",storm%track(1,1)
                print *,"   ERROR = ",i
                stop
            else if (i > storm%num_casts + 2) then
                print *,"Invalid storm indexing, i > num_casts + 2..."
                print *,"This really should not happen, what have you done?"
                stop
            endif
        endif

        ! Interpolate in time for all parameters
        if (i == storm%num_casts + 1) then
            i = i - 1
            ! At last forecast, use last data for storm strength parameters and
            ! velocity, location uses last velocity and constant motion forward

            ! Convert coordinates temporarily to meters so that we can use
            ! the pre-calculated m/s velocities from before
            x = latlon2xy(storm%track(2:3,i),storm%track(2:3,i))
            x = x + (t - storm%track(1,i)) * storm%velocity(:,i)
            
            fn = [xy2latlon(x,storm%track(2:3,i)), &
                  storm%velocity(:,i), storm%max_wind_radius(i), &
                  storm%max_wind_speed(i), storm%central_pressure(i), &
                  storm%rrp(i)]
        else
            ! Inbetween two forecast time points (the function storm_index
            ! ensures that we are not before the first data point, i.e. i > 1)
            tn = storm%track(1,i)
            tnm = storm%track(1,i-1)
            weight = (t - tnm) / (tn - tnm)
            fn = [storm%track(2:3,i),storm%velocity(:,i), &
                  storm%max_wind_radius(i),storm%max_wind_speed(i), &
                  storm%central_pressure(i), storm%rrp(i)]
            fnm = [storm%track(2:3,i - 1),storm%velocity(:,i - 1), &
                   storm%max_wind_radius(i - 1),storm%max_wind_speed(i - 1), &
                  storm%central_pressure(i - 1), storm%rrp(i - 1)]
            fn = weight * (fn - fnm) + fnm
        endif

        ! Set output variables
        location = fn(1:2)
        velocity = fn(3:4)
        max_wind_radius = fn(5)
        max_wind_speed = fn(6)
        central_pressure = fn(7)
        rrp = fn(8)

    end subroutine get_cle_storm_data
   
 
    ! ==========================================================================
    !  inner_derivative(f,r_m,v_m,r_a) 
    !   Evaluate the derivative of the inner model with respect to the 
    !   merge point radius. Fortran can not handle the computation of the 
    !   derivative in one line, thus it is calculated by breaking it up. 
    !   The analytical solution of the derivative is very messy, but possible 
    !   to solve.  
    ! ==========================================================================
    real(kind=8) function evaluate_inner_derivative(f,r_m,v_m,r_a,v_a) result(dMa)
        
        ! Input 
        real(kind=8), intent(in) :: f, r_m, v_m, r_a, v_a 
        
        ! Variables 
        real(kind=8) :: denominator, M_a  

        ! Calculate necessary components of the derivative to make 
        ! calculations easier  
        denominator = 2 - alpha + alpha * (r_a / r_m)**2
        M_a = 0.5 * f * r_a**2 + r_a * v_a 
        dMa = (2 * M_a / r_a) * (1/denominator) 
     
    end function evaluate_inner_derivative 

    ! ==========================================================================
    !  evaluate_v_a(f,r_m,v_m,r_a)  
    !    Evaluates the inner velocity at the merge point that allows the inner
    !    and outer model to be continuous. 
    ! ==========================================================================
    real(kind=8) function evaluate_v_a(f,r_m,v_m,r_a) result(v_a)
        
        ! Input 
        real(kind=8), intent(in) :: f, r_m, v_m, r_a

        v_a = ((2.0*(r_a/r_m)**2)/(2.0-alpha+alpha*(r_a/r_m)**2))**(1.0/(2.0-alpha))
        v_a = v_a*(0.5*f*r_m**2 + r_m*v_m) 
        v_a = v_a - 0.5*f*r_a**2
        v_a = v_a/r_a  
        
    end function evaluate_v_a  


    ! ==========================================================================
    !  solve_r_0(f,r_a,v_a,res,r_0, r_guess, m_out)  
    !    Determines the wind profile for a hurricane using 
    !    the inner and outer model for the Chavas model. 
    ! ==========================================================================
    subroutine solve_r_0(f, r_a, v_a, res, r_0, r_guess)
    
        ! Input 
        real(kind=8), intent(in) :: f, r_a, v_a
        integer, intent(in) :: res
        
        ! Variables with less strict i/o  
        real(kind=8) :: r_0, r_guess, m_diff, r_step, dr,r_m 
        real(kind=8) :: m_out(1:res)
        integer :: i  

        ! Parameters  
        real(kind=8), PARAMETER :: R_TOL=0.1  
        integer :: STEP_FLAG=1 

        ! Intialize the guess        
        r_0 = r_guess 
        r_step = r_guess

        ! Golden section method to determine r_0, the radius 
        ! at which the wind's speed drops to zero.
        do 
            call integrate_m_out(f, r_0, r_a, res, m_out)
            if (m_out(1) - 0.5*f*r_a**2 - r_a*v_a < 0) then
                r_step = r_step/STEP_FLAG  
                r_0 = r_0 + r_step   
            else
                r_step = r_step/2.0 
                r_0 = r_0 - r_step 
                STEP_FLAG = 2 
            end if
            
            ! Check the tolerance value 
            if (r_step < R_TOL) then 
                r_guess = r_0
                exit
            end if
        end do  
        
    end subroutine solve_r_0     

    ! ==========================================================================
    !  integrate_m_out(f, r_0, r_a, res, m_out) 
    !   Integrates from r_0 to r_a to determine the outer wind profile. 
    ! ==========================================================================
    subroutine integrate_m_out(f,r_0,r_a,res,m_out) 

        ! Input variables 
        real(kind=8), intent(in) :: f, r_0, r_a
        real(kind=8), dimension(1:res) :: m_out
        integer, intent(in) :: res 

        ! Parameters and Other variables 
        real(kind=8) :: V_m, dr, r_guess, r_p 
        integer :: i  

        ! Intialize 
        m_out(res) = 0.5*f*r_0**2 
        dr = (r_0 - r_a)/res  
        r_p = r_0 - dr 
        m_out(res-1) = 0.5*f*r_p**2 + r_p*f*(r_0-r_p)
 
        ! Integrates the outer model equation from r_0 to r_a 
        do i = res-2, 1, -1
            r_guess = r_p
            r_p = r_guess - dr 
            V_m = (m_out(i+1) - 0.5*f*r_p**2)/r_p 
            m_out(i) = m_out(i+1) - chi*dr*((r_p*V_m)**2)/(r_0**2 - r_p**2) 
        end do 
 
    end subroutine integrate_m_out

 
    ! ==========================================================================
    !  solve_hurricane_wind_parameters(f, r_m, v_m, res, r_0, r_a)
    !   Determine r_0, r_a, and v_a.  
    !   Determine the merge point when given the f value. At this current 
    !   moment we are using values of chi = 1.0 and alpha = 1.0. 
    ! ==========================================================================
    subroutine solve_hurricane_wind_parameters(f, r_m, v_m, res, r_0, r_a)
    
        ! Input 
        real(kind=8), intent(in) :: f, r_m, v_m
        integer, intent(in) :: res 
    
        ! Variables 
        real(kind=8) :: r_0, r_0_guess, r_a, r_step, v_a  
        real(kind=8) :: slope_difference
        real(kind=8) :: m_out(1:res) 
        
        ! Parameters
        real(kind=8), parameter :: r_tol = 0.1  
        integer :: step_flag = 1, i 

        ! Parameters
        real(kind=8) :: dr, r
        real(kind=8) :: inner_res, outer_res 
        real(kind=8), dimension(1:res) :: v_out, v_in  
        
        ! Initialize guesses for merge and r_0 
        r_a = 2.0*r_m 
        r_0_guess = 5.0*r_m  
        r_step = r_a  
    
        do 
            v_a = evaluate_v_a(f,r_m,v_m,r_a)
            !print *, "v_a: ", v_a 
            if (v_a < 0) then 
                r_step = r_step/2 
                r_a = r_a - r_step 
                cycle 
            end if 
            
            call solve_r_0(f,r_a,v_a,res,r_0,r_0_guess)  
            
            slope_difference = evaluate_inner_derivative(f,r_m,v_m,r_a,v_a) - & 
                                chi*((r_a*v_a)**2)/(r_0**2 - r_a**2) 
            if (slope_difference < 0) then 
                r_step = r_step/2.0 
                r_a = r_a - r_step
                step_flag = 1  
            else
                r_step = r_step/step_flag 
                r_a = r_a + r_step 
            end if  
            
            if (r_step < r_tol) then 
                exit 
            end if 
        end do  
 
    end subroutine solve_hurricane_wind_parameters  

    ! ==========================================================================
    !  solve_wind_profile(f, r_m, v_m, res) 
    ! ==========================================================================
    subroutine solve_wind_profile(f, r_m, v_m, res)
        
        ! Input 
        real(kind=8), intent(in) :: f, r_m, v_m
        integer, intent(in) :: res 
        real(kind=8), save :: r_0, r_a, v_a 

        ! Parameters
        real(kind=8) :: dr, r
        real(kind=8) :: inner_res, outer_res 
        real(kind=8), dimension(1:res) :: v_out, v_in, m_out  
        integer :: i  

        call solve_hurricane_wind_parameters(f, r_m, v_m, res, r_0, r_a)
        
        ! Determine the inner model velocity vector 
        dr = (r_a/res)
        r = dr  
        v_in(1) = 0 
        
  
        do i = 2, res, 1
            v_in(i) = evaluate_v_a(f, r_m, v_m, r) 
            r = r + dr
        end do
        
        ! Determine the outer model velocity vector 
        dr = (r_0 - r_a)/res
        r = r_0 - dr 
        v_out(res) = 0
        
        call integrate_m_out(f, r_0, r_a, res, m_out) 
 
        do i = res-1, 1, -1 
            v_out = (m_out(i) - 0.5*f*r**2)/r 
            r = r - dr 
        end do 
        
    end subroutine solve_wind_profile 


    ! ==========================================================================
    !  set_cle_storm_fields()
    ! ==========================================================================
    subroutine set_cle_storm_fields(maux,mbc,mx,my,xlower, &
                                    ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, storm)

        use geoclaw_module, only: g => grav, ambient_pressure
        use geoclaw_module, only: coriolis, deg2rad
        use geoclaw_module, only: spherical_distance

        use geoclaw_module, only: rad2deg

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(cle_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2) 
        real(kind=8) :: f, mwr, mws, Pc, Pa, dp, wind, tv(2), rrp 
        real(kind=8) :: mod_mws, trans_speed, ramp 
        integer :: i,j 
        
        ! Variables for CLE model calculations
        integer, PARAMETER :: res=2000
        real(kind=8), save :: v_vec(1:res), p_vec(1:res)
        real(kind=8), save :: m_out_vec(1:res)
        real(kind=8), save :: last_time=-1
        real(kind=8), save :: dr, r_a, r_0
        integer :: io_status
        integer, save :: in_res, out_res

        
        ! Additional local storage for CLE Model
        real(kind=8) :: v, next_loc(2)
        integer :: r_index
        real(kind=8) :: r_factor, t_factor

        ! Output the pressure and velocity fields.  
        ! Determines the radial wind ONLY IF called at a new time. 

        !First, interpolate all of the relevant tracked storm parameters.
        call get_cle_storm_data(t,storm,sloc,tv,mwr,mws,Pc,rrp)

        ! Other quantities of interest 
        Pa = ambient_pressure 

        ! Additional quantities of interest 
        i = storm_index(t, storm)
        f = coriolis(sloc(2))
        !tv = storm%velocity(:, i)
        !Pa = ambient_pressure
        !sloc = cle_storm_location(t, storm)
        !f = coriolis(sloc(2))
        !Pc = storm%central_pressure(i)
        !mws = storm%max_wind_speed(i)
        !mwr = storm%max_wind_radius(i)

        ! Calculate CLE parameters 
        ! Subtract translational speed of storm from maximum wind speed 
        ! to avoid distorition in the CLE curve fit. Added back later
        trans_speed = sqrt(tv(1)**2 + tv(2)**2) 
        mod_mws = mws - trans_speed 
    
        ! Convert wind speed (10 m) to top of atmospheric boundary layer 
        mod_mws = mws - trans_speed 
        
        ! Fill array with CLE model parameters 
        if (last_time /= t) then
            ! Determine the wind profile, but only if the profile isn't already
            ! saved. In v_vec and p_vec. 
            last_time = t
            call solve_hurricane_wind_parameters(f,mwr,mws,res,r_0,r_a)
            dr = real(r_0) / (res - 1)
            !out_res = number of points in outer model/in_res = points in inner model
            out_res = ceiling(real((res - 1) * (r_0 - r_a)) / r_0)
            in_res = res - out_res
            !Determine angular momentum profile for outer model
            call integrate_m_out(f,r_0,r_0 - out_res*dr,out_res,m_out_vec)
            do j= 1, out_res
               !convert angular momentum to azimuthal velocity.
               r = r_0 - (j-1)*r_0/res
               m_out_vec(out_res + 1 - j) = (m_out_vec(out_res + 1 - j)- .5*f*r**2)/r
            end do

            p_vec(1) = 0
            v_vec(1) = 0
           
            do j=2, res
                ! Combine inner and outer model into single radial wind profile.
                if(j <= in_res) then
                   r = (j - 1)*dr
                   v_vec(j) = evaluate_v_a(f, mwr, mws, r)
                else
                   v_vec(j) = m_out_vec(j - in_res)
                end if
                ! Pressure is determined using the gradient wind equation.
                p_vec(j) = p_vec(j - 1) + (v_vec(j)**2)/r + f*v_vec(j)
            end do
           
            do j=1, res
                ! Normalize pressure to match measurements.
                p_vec(j) = (Pa - Pc)*p_vec(j)/p_vec(res) + Pc
            end do
            
            open(unit=4242, file='test.txt',status='replace', & 
                 action='write',iostat=io_status) 
            write(4242,*) v_vec 
        end if
        
        ! Set fields 
        do j=1-mbc,my+mbc
           ! Assign a clockwise, radially symmetric wind field by piece-wise linear
           ! interpolation of v_vec.  Likewise with p_vec.
            y = ylower + (j-0.5d0) * dy ! Degrees latitude
            f = coriolis(y)
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx ! Degrees longitude 

                ! Calcuate storm centric polar coordinate location of grid
                ! cell center, uses Haversine formula
                r = spherical_distance(x, y, sloc(1), sloc(2))
                r_index = floor(r/dr)
                r_factor = (r - r_index * dr) / dr
                r_index = r_index + 1 
                theta = atan2((y - sloc(2)) * DEG2RAD,(x-sloc(1))*DEG2RAD)
                if (r_index < res) then
                    aux(pressure_index, i, j) = (1 - r_factor) * p_vec(r_index) & 
                                        + r_factor * p_vec(r_index + 1) 
                    aux(wind_index, i, j) = -0.9 * ((1 - r_factor) * v_vec(r_index) & 
                                        + r_factor * v_vec(r_index + 1)) * sin(theta)  
                    aux(wind_index+1, i, j) = 0.9 * ((1 - r_factor) * v_vec(r_index) & 
                                        + r_factor * v_vec(r_index + 1)) * cos(theta)  
                else
                    aux(pressure_index,i,j) = Pa
                    aux(wind_index,i,j) = 0
                    aux(wind_index+1,i,j) = 0
                end if
            end do
        end do

    end subroutine set_cle_storm_fields
 
end module cle_storm_module

!!!!!! - Deprecated code !!!!!

!    ! ==========================================================================
!    !  set_cle_storm_fields()
!    ! ==========================================================================
!    subroutine set_cle_storm_fields(maux,mbc,mx,my,xlower, &
!                                    ylower,dx,dy,t,aux, wind_index, &
!                                    pressure_index, storm)
!
!        use geoclaw_module, only: g => grav, ambient_pressure
!        use geoclaw_module, only: coriolis, deg2rad
!        use geoclaw_module, only: spherical_distance
!
!        use geoclaw_module, only: rad2deg
!
!        implicit none
!
!        ! Time of the wind field requested
!        integer, intent(in) :: maux,mbc,mx,my
!        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t
!
!        ! Storm description, need in out here since we may update the storm
!        ! if at next time point
!        type(cle_storm_type), intent(in out) :: storm
!
!        ! Array storing wind and presure field
!        integer, intent(in) :: wind_index, pressure_index
!        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
!        integer, PARAMETER :: res=2000
!        real(kind=8), save :: v_vec(1:res), p_vec(1:res)
!        real(kind=8), save :: m_out_vec(1:res)
!        real(kind=8), save :: last_time=-1
!        real(kind=8), save :: dr, r_a, r_0
!        integer :: io_status
!        integer, save :: in_res, out_res
!
!        ! Local storage
!        real(kind=8) :: x, y, r, theta, sloc(2), v, next_loc(2)
!        real(kind=8) :: f, mwr, mws, Pc, Pa, dp, wind, tv(2)
!        integer :: i,j, r_index
!        real(kind=8) :: r_factor, t_factor
!
!        !output the pressure and velocity fields.  Determines the radial wind
!        !profile only if the called at a new time.
!
!        !First, interpolate all of the relevant tracked storm parameters.
!        i = storm_index(t, storm)
!        tv = storm%velocity(:, i)
!        Pa = ambient_pressure
!        sloc = cle_storm_location(t, storm)
!        f = coriolis(sloc(2))
!        Pc = storm%central_pressure(i)
!        mws = storm%max_wind_speed(i)
!        mwr = storm%max_wind_radius(i)
!        
!        if (last_time /= t) then
!           !Determine the wind profile, but only if the profile isn't already
!           !saved. In v_vec and p_vec
!           last_time = t
!           call solve_hurricane_wind_parameters(f, mwr, mws, res, r_0, r_a)
!           dr = real(r_0)/(res - 1)
!           !out_res = number of points in outer model/in_res = points in inner model
!           out_res = ceiling( real((res - 1)*(r_0 - r_a))/r_0 )
!           in_res = res - out_res
!           !Determine angular momentum profile for outer model
!           call integrate_m_out(f, r_0, r_0 - out_res*dr, out_res, m_out_vec)
!           do j= 1, out_res
!              !convert angular momentum to azimuthal velocity.
!              r = r_0 - (j-1)*r_0/res
!              m_out_vec(out_res + 1 - j) = (m_out_vec(out_res + 1 - j)- .5*f*r**2)/r
!
!           end do
!
!           p_vec(1) = 0
!           v_vec(1) = 0
!           do j=2, res
!              !combine inner and outer model into single radial wind profile.
!              if(j <= in_res) then
!                 r = (j - 1)*dr
!                 v_vec(j) = evaluate_v_a(f, mwr, mws, r)
!              else
!                 v_vec(j) = m_out_vec(j - in_res)
!              end if
!              !pressure is determined using the gradient wind equation.
!              p_vec(j) = p_vec(j - 1) + (v_vec(j)**2)/r + f*v_vec(j)
!           end do
!           do j=1, res
!              !normalize pressure to match measurements.
!              p_vec(j) = (Pa - Pc)*p_vec(j)/p_vec(res) + Pc
!           end do
!            
!           open(unit=4242, file='test.txt',status='replace', & 
!                action='write',iostat=io_status) 
!           write(4242,*) v_vec 
!        end if
!        
!        ! Set fields 
!        do j=1-mbc,my+mbc
!           !Assigns a clockwise, radially symmetric wind field by piece-wise linear
!           !interpolation of v_vec.  Likewise with p_vec.
!            y = ylower + (j-0.5d0) * dy ! Degrees latitude
!            f = coriolis(y)
!            do i=1-mbc,mx+mbc
!                x = xlower + (i-0.5d0) * dx  
!
! 
!            r = spherical_distance(x, y, sloc(1), sloc(2))
!            r_index = floor(r/dr)
!            r_factor = (r - r_index*dr)/dr
!            r_index = r_index + 1
!            theta = atan2((y - sloc(2)) * DEG2RAD,(x - sloc(1)) * DEG2RAD)
!            if(r_index < res) then
!              aux(pressure_index, i, j) = (1-r_factor)*p_vec(r_index) &
!                                        + r_factor*p_vec(r_index + 1)
!              aux(wind_index,i,j) = -0.9*( (1-r_factor)*v_vec(r_index) + &
!                                       r_factor*v_vec(r_index + 1) )*sin(theta)
!              aux(wind_index+1,i,j) = 0.9*( (1-r_factor)*v_vec(r_index)  +&
!                                       r_factor*v_vec(r_index + 1) )*cos(theta)
!            else
!              aux(pressure_index,i,j) = Pa
!              aux(wind_index,i,j) = 0
!              aux(wind_index+1,i,j) = 0
!            end if
!          end do
!        end do
!
!    end subroutine set_cle_storm_fields
!
! 
!end module cle_storm_module
!
!
