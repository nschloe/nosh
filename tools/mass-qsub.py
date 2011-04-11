#!/usr/bin/env python
'''Mass-submits jobs for timings on machines using qsub. Intended for
timing jobs.'''
# ==============================================================================
import sys, os
# ==============================================================================
def _main():

    import datetime, math, socket

    numprocs_range = [1]
    numprocs_range.extend( range( 8, 257, 8 ) )

    queue       = "qshort"
    walltime    = "0:05:00"
    executable  = "$HOME/code/ginla/build/mpi/scalability-tests/scaltest.exe"
    options     = "" #"--input=init/cutcircle1000-balanced.par"
    name        = "epetra-generic" # "belos-timer"
   
    # --------------------------------------------------------------------------
    # Create the submit scripts.
    # This could possibly be done with a command line call too, but
    # this way it is possible to check up on the parameters.
    hostname = socket.gethostbyaddr(socket.gethostname())[0]
    date = datetime.datetime.now().isoformat(' ')
    submit_files = []

    # Use either
    # harpertown:ib (that's 4 cores per CPU times 2 for hyperthreading)
    # or
    # westmere (that's 12 cores per CPU times 2 for hyperthreading).
    # All westmere nodes are InfiniBand (ib) anyways.
    # node_type = "harpertown:ib"
    # num_cores_per_node = 8
    # Prefer Westmere b/c as of now (Jan 25 2011) they represent the largest
    # homogenous set of processes on CalcUA.
    node_type = "westmere"
    num_cores_per_node = 24
    for numprocs in numprocs_range:

        # Compute the number of nodes needed to fir numprocs processes.
        if numprocs % num_cores_per_node == 0:
            num_nodes = numprocs / num_cores_per_node
        else:
            num_nodes = math.ceil( numprocs / float(num_cores_per_node) )

        jobname = "%s-core%03d-nodes%03d-ppn%02d" % ( name, numprocs, num_nodes, num_cores_per_node )

        # Open the submit file.
        submitfile_name = "%s.pbs" % jobname
        submit_files.append( submitfile_name )

        file_handle = open( submitfile_name, 'w' )

        stdout_file = "out/%s.o$PBS_JOBID" % jobname
        stderr_file = "out/%s.e$PBS_JOBID" % jobname

        resource_list = "nodes=%d:ppn=%d:%s" \
                        % ( num_nodes, num_cores_per_node, node_type )

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
#PBS -l %s

  # make sure we are in the right directory in case writing files 
cd $PBS_O_WORKDIR

  # load the environment
module load ictce
module load scripts

# Prepend some info to the output file
echo "# host: %s"
echo "# arch: %s"
echo "# numprocs: %d"

mympirun %s %s''' \
% ( hostname, date,
    jobname,
    queue,
    walltime,
    stdout_file,
    stderr_file,
    resource_list,
    hostname,
    resource_list,
    numprocs,
    executable, options
  )
)

        # Close file. 
        file_handle.close()
    # --------------------------------------------------------------------------
    # submit the jobs
    num_samples = 10
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
