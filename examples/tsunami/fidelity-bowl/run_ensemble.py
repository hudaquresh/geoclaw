#!/usr/bin/env python

"""Run specified faults in the provided CSV file."""

from __future__ import print_function

import sys
import os

import numpy
import matplotlib.pyplot as plt

import batch.habanero as batch

import clawpack.geoclaw.dtopotools as dtopotools


class FidelityJob(batch.HabaneroJob):

    r"""Job describing a single Okada based fault relization.

    Fault parameterization is:

    """

    def __init__(self):
        r"""
        Initialize a FidelityJob object.

        See :class:`FidelityJob` for full documentation

        """

        super(FidelityJob, self).__init__()

        self.type = "bowl-radial"
        self.name = "fidelity-test"
        self.prefix = "fault_%s" % self.run_number
        self.executable = 'xgeoclaw'

        # Data objects
        import setrun
        self.rundata = setrun.setrun()

        # No variable friction for the time being
        self.rundata.friction_data.variable_friction = False

        # Replace dtopo file with our own
        self.dtopo_path = 'fault_s%s.tt3' % self.base_subfault.slip
        self.rundata.dtopo_data.dtopofiles = [[3, 4, 4, self.dtopo_path]]


    def __str__(self):
        output = super(FaultJob, self).__str__()
        output += "\n  Fault Parameters:\n"
        output += "      run_number = %s\n" % self.run_number
        output += "      slips = "
        for subfault in self.fault.subfaults:
            output += "%s " % subfault.slip
        output += "\n"
        return output


    def write_data_objects(self):

        # Create dtopo file
        x, y = self.fault.create_dtopo_xy()
        dtopo = self.fault.create_dtopography(x, y)
        dtopo.write(path=self.dtopo_path, dtopo_type=3)

        # Plot fault here
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)

        cmap = plt.get_cmap("YlOrRd")
        self.fault.plot_subfaults(axes=axes, slip_color=True, cmap_slip=cmap, 
                                  cmin_slip=FaultJob.cmin_slip,
cmax_slip=FaultJob.cmax_slip,
                                  plot_rake=True)
        axes.set_title("$M_o = %s$, $M_w = %s$" % (str(self.fault.Mo()),
str(self.fault.Mw())))
        fig.savefig("fault_slip.png")

        # Write other data files
        super(FaultJob, self).write_data_objects()


if __name__ == '__main__':

    # If a file is found on the command line, use that to calculate the 
    # quadrature points, otherwise calculate it given a default range
    if len(sys.argv) > 1:
        path = sys.argv[1]

        # Load fault parameters
        slips = numpy.loadtxt(path)
        if len(slips.shape) == 1:
            slips = [slips]
    else:
        slips = numpy.loadtxt("./random_sample.txt")
        # slip_range = (0.0, 120.0)
        # slips = calculate_quadrature(slip_range)
    
    # Create all jobs
    path = os.path.join(os.environ.get('DATA_PATH', os.getcwd()), 
                "tohoku", "okada-fault-random",
                    "run_log.txt")
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path))

    # Find minimum and maximum slips for plotting of the faults
    FaultJob.cmin_slip = numpy.min(slips)
    FaultJob.cmax_slip = numpy.max(slips)

    with open(path, 'w') as run_log_file: 
        jobs = []
        for (n, slip) in enumerate(slips):
            run_log_file.write("%s %s\n" % (n, ' '.join([str(x) for x in
slip])))
            jobs.append(FaultJob(slip, run_number=n))

    controller = batch.HabaneroBatchController(jobs)
    controller.wait = False
    controller.plot = False
    print(controller)
    # controller.run()
