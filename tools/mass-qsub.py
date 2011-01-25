#!/usr/bin/env python
'''Mass-submits jobs for timings on machines using qsub. Intended for
timing jobs.'''
# ==============================================================================
import sys, os
# ==============================================================================
def _main():

    import datetime

    numprocs_range = range( 8, 257, 8 )

    queue       = "qshort"
    walltime    = "0:05:00"
    executable  = "$HOME/code/ginla/build/mpi/packages/ginla-fvm/examples/linear-solve-test/keo-matvec.exe"
    options     = "--input=init/cutcircle1000-balanced.par"
   
    # --------------------------------------------------------------------------
    # Create the submit scripts.
    # This could possibly be done with a command line call too, but
    # this way it is possible to check up on the parameters.
    hostname = os.uname()[1]
    date = datetime.datetime.now().isoformat(' ')
    submit_files = []
    for numprocs in numprocs_range:
        num_nodes          = numprocs / 8
        num_cores_per_node = 8

        jobname = "belos-timer-core%03d" % numprocs

        # Open the submit file.
        submitfile_name = "%s.pbs" % jobname
        submit_files.append( submitfile_name )

        file_handle = open( submitfile_name, 'w' )

        stdout_file = "out/%s.o$PBS_JOBID" % jobname
        stderr_file = "out/%s.e$PBS_JOBID" % jobname

        # Write contents.
        file_handle.write( '''#/bin/bash
# ------------------------------------------------------------------------------
# submit file for %s, %s
# ------------------------------------------------------------------------------

#PBS -N %s
#PBS -q %s
#PBS -l walltime=%s

#PBS -o %s
#PBS -e %s

  # Specify number of cores.
  # Use Westmere here b/c as of now (Jan 25 2011) they represent the largest
  # homogenous set of processes on CalcUA.
#PBS -l nodes=%d:ppn=%d:ib:westmere

  # make sure we are in the right directory in case writing files 
cd $PBS_O_WORKDIR

  # load the environment
module load ictce
module load scripts

mympirun %s %s''' \
% ( hostname, date,
    jobname,
    queue,
    walltime,
    stdout_file,
    stderr_file,
    num_nodes, num_cores_per_node,
    executable, options
  )
)

        # Close file. 
        file_handle.close()
    # --------------------------------------------------------------------------
    # submit the jobs
    num_samples = 1
    for submit_file in submit_files:
        for k in xrange( num_samples ):
            output = _run( "qsub %s" % submit_file )
    # --------------------------------------------------------------------------
    return 0
# ==============================================================================
def _parse_options():
    '''Parse input options.'''
    import optparse

    usage = "usage: %prog filename"

    parser = optparse.OptionParser( usage = usage )

    (options, args) = parser.parse_args()

    return args[0]
# ==============================================================================
def _run( command ):
    """Runs a given command on the command line and returns its output.
    """
    import subprocess

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
        sys.exit( "\nERROR: The command \n\n%s\n\nreturned a nonzero " \
                  "exit status. The error message is \n\n%s\n\n" \
                  "Abort.\n" % \
                  ( command, process.stderr.read()[:-1] )
                )
    return output
# ==============================================================================
if __name__ == '__main__':
    status = _main()
    sys.exit( status )
# ==============================================================================
