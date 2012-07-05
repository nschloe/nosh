#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Function to plot the continuation curves with stability information."""
# ==============================================================================
import subprocess
import os
import matplotlib.pyplot as pp
import matplotlib.image as mpimg
import numpy as np
# ==============================================================================
def _main():

    args = _parse_args()
    continuation_data = np.loadtxt( args.continuation_data )

    # load image
    for k in xrange(len(args.stateimgs)):
        fig = pp.figure(figsize=(12.0, 5.0))
        ax  = fig.add_subplot(1, 2, 1)
        # Plot parameter data.
        plot_columns = [3, 4]
        x_values = continuation_data[:,plot_columns[0]]
        y_values = continuation_data[:,plot_columns[1]]
        pp.plot( x_values, y_values, 'k-' )
        pp.xlabel('mu')
        pp.ylabel('Gibbs energy')
        pp.ylim( -1.0, -0.0 )
        pp.plot(x_values[k], y_values[k], 'or')

        # trim the image
        _run('convert -trim %s tmp.png' % args.stateimgs[k])
        ax  = fig.add_subplot(1, 2, 2)
        img = mpimg.imread('tmp.png')
        imgplot = pp.imshow(img)
        # hide axis
        ax.set_frame_on(False)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        outfile = 'diptych%04d.png' % k
        pp.savefig(outfile)
        #pp.show()
    # remove temp file
    os.remove('tmp.png')

    # create movie
    _run('avconv -r 5 -i diptych%04d.png -f webm test.webm')

    return
# ==============================================================================
def _run( command ):
    """Runs a given command on the command line and returns its output.
    """
    print command

    process = subprocess.Popen( command,
                                shell     = True,
                                stdout    = subprocess.PIPE,
                                stderr    = subprocess.PIPE,
                                close_fds = True
                              )
    output = process.stdout.read()[:-1]
    ret = process.wait()

    if ret != 0:
        import sys
        sys.exit( "\nERROR: The command \n\n%s\n\nreturned a nonzero " \
                  "exit status. The error message is \n\n%s\n\n" \
                  "Abort.\n" % \
                  ( command, process.stderr.read()[:-1] )
                )
    return output
# ==============================================================================
def _parse_args():
    '''Parse input arguments.'''
    import argparse

    parser = argparse.ArgumentParser( description = 'Create animated plots.' )

    parser.add_argument( 'stateimgs',
                         metavar = 'FILES',
                         nargs = '+',
                         help = 'images from states'
                       )

    parser.add_argument( '--continuation-data', '-c',
                         required = True,
                         type = str,
                         help = 'continuation data file'
                       )

#    parser.add_argument('--timesteps', '-t',
#                         metavar='TIMESTEP',
#                         type=int,
#                         nargs = '+',
#                         default=None,
#                         help='read a particular time step/range (default: all)'
#                        )

    return parser.parse_args()
# ==============================================================================
if __name__ == '__main__':
    _main()
# ==============================================================================
