
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw.topotools import topo1writer, topo2writer
from clawpack.geoclaw import dtopotools
from numpy import *
import os


def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints=201
    nypoints=201
    xlower = -10e0
    xupper= 10e0
    ylower = -10e0
    yupper= 10e0
    outfile= "topo1.topotype2"
    topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def maketopo2():
    """
    Another topo file for testing
    """
    nxpoints=21
    nypoints=51
    xlower = -0.02
    xupper= 0.18
    ylower = -0.1
    yupper= 0.4
    outfile= "topo2.topotype2"
    topo2writer(outfile,topo2,xlower,xupper,ylower,yupper,nxpoints,nypoints)


def topo(x,y):
    z = -10*ones(x.shape)
    return z

def topo2(x,y):
    z = -15*ones(x.shape)
    return z


def read_subfaults_dtopo1(plotfig=None):
    """
    Test data
    """
    fname_subfaults = 'dtopo1.csv'

    # Format of subfault file:
    columns = """longitude latitude depth length width strike dip rake slip
       """.split()
    defaults = {'latlong_location': 'top center'}
    units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}

    subfaults = dtopotools.read_subfault_model(fname_subfaults, \
                        columns=columns, units=units, \
                        defaults = defaults, skiprows=1, delimiter=',')
                    
    if plotfig:
        figure(plotfig)
        dtopotools.plot_subfaults(subfaults,slip_color=True,plot_rake=True)
        
    return subfaults

def make_dtopo_dtopo1(plotfig=None):
    """
    Test data.
    """
    from clawpack.geoclaw import okada2

    subfaults = read_subfaults_dtopo1()

    dtopo_fname = 'dtopo1.tt3'

    if os.path.exists(dtopo_fname):
        print "*** Not regenerating dtopo file (already exists): %s" \
                % dtopo_fname
    else:
        print "Using Okada model to create %s " % dtopo_fname

        # Needed for extent of dtopo file:
        xlower = -0.4
        xupper = 0.6
        ylower = -0.4
        yupper = 0.4

        # number of grid points in dtopo file
        mx = 51
        my = 41

        # Create dtopo_params dictionary with parameters for dtopo file: 
        dtopo_params = {}
        dtopo_params['fname'] = dtopo_fname
        dtopo_params['faulttype'] = 'static'
        dtopo_params['dtopotype'] = 3
        dtopo_params['mx'] = mx
        dtopo_params['my'] = my
        dtopo_params['xlower'] = xlower
        dtopo_params['xupper'] = xupper
        dtopo_params['ylower'] = ylower
        dtopo_params['yupper'] = yupper
        dtopo_params['t0'] = 0.
        dtopo_params['tfinal'] = 100.
        dtopo_params['ntimes'] = 5

        dtopo = dtopotools.make_dtopo_from_subfaults(subfaults, dtopo_params)

        if plotfig:
            figure(plotfig)
            x = dtopo.x
            y = dtopo.y
            dz_final = dtopo.dz_list[-1]
            dtopotools.plot_dz_colors(x,y,dz_final,cmax_dz=16,dz_interval=1)

        return dtopo

if __name__=='__main__':
    maketopo()
    maketopo2()
    make_dtopo_dtopo1()
