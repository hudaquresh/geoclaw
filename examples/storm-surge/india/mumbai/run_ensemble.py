#!/usr/bin/env python

"""Run specified faults in the provided CSV file."""

from __future__ import print_function

import sys
import os

import numpy
import matplotlib.pyplot as plt

import batch.habanero as batch

import clawpack.geoclaw.dtopotools as dtopotools
import datetime
import time 

 
def days2seconds(days):
    return days * 60.0**2 * 24.0

class MumbaiJob(batch.HabaneroJob):

    r"""Job describing a Mumbai run .

    """

    def __init__(self, run_number, t0, tf, storm_dir):  
        r"""
        Initialize a MumbaiJob object.

        See :class:`MumbaiJob` for full documentation

        """

        super(MumbaiJob, self).__init__()

        self.type = "storm-surge"
        self.name = "h80-amr5"
        self.run_number = run_number
        self.storm_dir = storm_dir 
        self.t0 = datetime.datetime(t0[0],t0[1],t0[2],t0[3]) - \
                        datetime.datetime(t0[0], 1, 1, 0)    
        self.tf = datetime.datetime(tf[0],tf[1],tf[2],tf[3]) - \
                        datetime.datetime(t0[0],t0[1],t0[2],t0[3])      
        self.prefix = "%s" % self.run_number
        self.executable = 'xgeoclaw'

        # Data objects
        import setrun
        self.rundata = setrun.setrun()

 
        # Set rundata friction, surge, and clawdata variables  
        self.rundata.friction_data.variable_friction = True

        self.rundata.surge_data.storm_type = 1      # H80
        #self.rundata.surge_data.storm_type = 4      # CLE 
        #self.rundata.surge_data.storm_type = 5      # H10 
        #self.rundata.surge_data.storm_type = 6      # SLOSH
 
        self.rundata.clawdata.t0 = days2seconds(self.t0.days) + \
                                                self.t0.seconds
 
        self.rundata.clawdata.tfinal = days2seconds(self.t0.days + \
                                                    self.tf.days)

        # == setregions.data values ==
        # Mumbai Region  
        self.rundata.regiondata.regions.append([2, 5, self.rundata.clawdata.t0, 
                                                      self.rundata.clawdata.tfinal,
                                                      70, 75, 17, 22])
        # Mumbai 
        self.rundata.regiondata.regions.append([4, 7, days2seconds(self.t0.days + 1.0) + self.t0.seconds, 
                                                      self.rundata.clawdata.tfinal,
                                                      72.6, 73, 18.80, 19.15])

        # Set Gauge Data 
        # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
        self.rundata.gaugedata.gauges = [] 
        self.rundata.gaugedata.gauges.append([1, 72.811790, 18.936508,
                                            self.rundata.clawdata.t0, 
                                            self.rundata.clawdata.tfinal])  
        self.rundata.gaugedata.gauges.append([2, 72.972316, 18.997762,
                                            self.rundata.clawdata.t0, 
                                            self.rundata.clawdata.tfinal])  
        self.rundata.gaugedata.gauges.append([3, 72.819311, 18.818044,
                                            self.rundata.clawdata.t0, 
                                            self.rundata.clawdata.tfinal])  
        self.rundata.gaugedata.gauges.append([4, 72.50, 18.50,
                                            self.rundata.clawdata.t0, 
                                            self.rundata.clawdata.tfinal]) 

        self.rundata.gaugedata.aux_out_fields = [4, 5, 6] 
 
                                                  
        # Set surge data  
        self.rundata.surge_data.display_landfall_time = True 
        self.rundata.surge_data.landfall = days2seconds(self.t0.days) + \
                                                        self.t0.seconds 
        storm_number = str(self.run_number)
        storm_file_no = 'mumbai_%s.storm' % storm_number
        self.rundata.surge_data.storm_file = os.path.join(self.storm_dir, storm_file_no)


    def __str__(self):
        output = super(MumbaiJob, self).__str__()
        output += "\n  Mumbai Storm:\n"
        output += "      run_number = %s\n" % self.run_number
        output += "\n"
        return output


def run_mumbai_job(storm_0=0, storm_end=2): 
    r"""
    Setup jobs to run at specific storm start and end.  
    """
    path = os.path.join(os.environ.get('DATA_PATH', os.getcwd()), 
                        "india", "mumbai-storm-%i-%i" %(storm_0, storm_end),
                        "run_log.txt")

    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path))

    mumbai_t0 = numpy.loadtxt("storms/mumbai.t0", dtype=int, delimiter=",",\
                                skiprows = 0)  
    mumbai_tf = numpy.loadtxt("storms/mumbai.tf", dtype=int, delimiter=",",\
                                skiprows = 0)  
    storms = os.path.join(os.getcwd(), "storms", "mumbai-tracks")  
 
    with open(path, 'w') as run_log_file: 
        jobs = []
        for n in range(storm_0, storm_end + 1):  
            storm_name = 'mumbai_%s.storm' % str(n)  
            run_log_file.write("%s %s\n" % (n, '%s' %storm_name))
            t_final = mumbai_tf[n]
            t_0 = mumbai_t0[n]
            jobs.append(MumbaiJob(run_number=n,t0=t_0,tf=t_final,storm_dir=storms)) 
    
            if n == storm_end:  
                controller = batch.HabaneroBatchController(jobs)
                controller.email = 'hq2152@columbia.edu'
                print(controller)
                controller.run()
                jobs = []
                break 


if __name__ == '__main__':
    r"""
    TODO: Batch pause jobs.   
    """
    
    # Default batch run 
    # run_mumbai_job() 
    
    ## Batch 1 (0-599) 
    #start = 0 
    #end = 599 
    #run_mumbai_job(storm_0=start, storm_end=end) 
 
    ## Batch 2 (600-1199) 
    #start = 600  
    #end = 1199 
    #run_mumbai_job(storm_0=start, storm_end=end) 

    # Batch 3 (1200-1849) 
    start = 1200 
    end = 1849 
    run_mumbai_job(storm_0=start, storm_end=end) 

 
