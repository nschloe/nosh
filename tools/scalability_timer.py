#!/usr/bin/env python
'''Performs scalability runs, i.e., runs a given code through mpiexec -n X
with various X, parses the output and stores it in a file. Typically, one
would store timings then of course.'''
# ==============================================================================
import sys, subprocess, re
# ==============================================================================
def _main():
        filename = "scalingdata.dat"
        f = open( filename, "w" )

        # run over the number of procs
        min_numprocs = 2
        max_numprocs = 6

        num_runs = 5

        # write header
        header_string = "#\t"
        for num_procs in xrange( min_numprocs, max_numprocs+1 ):
            header_string += "%d\t\t" % num_procs
        header_string += "\n"
        f.write( header_string )


        for k in xrange( 0, num_runs ):
            for num_procs in xrange( min_numprocs, max_numprocs+1 ):
                max_time = _testrun( num_procs )
                f.write( "\t%e" % max_time )
            f.write( "\n" )

        f.close()

	return
# ==============================================================================
def _testrun( num_procs ):

        #mpiexec_cmd = "mpiexec"
        #test_exe = "/home/nico/ginla/build/mpi/packages/ginla-fvm/examples/linear-solve-test/keo-belos.exe"
        #test_options = "--input=cutcircle300-balanced.par"

        #cmd = [ mpiexec_cmd, "-n %d" % num_procs, test_exe, test_options ]
        cmd = "mpiexec -n %d /home/nico/ginla/build/mpi/packages/ginla-fvm/examples/linear-solve-test/keo-belos.exe --input=cutcircle300-balanced.par" % num_procs

        # run the test
        output = _run( cmd )

        # Parse the output, get the maximum execution time over all processes.
        # Find something of the kind
        # "Belos: PseudoBlockCGSolMgr total solve time    0.3433 (1)        0.3433 (1)        0.3434 (1)        "
        # Enclose the three timings at the in parentheses
        key = "Belos: PseudoBlockCGSolMgr total solve time"
        regex = "%s\s*(\d\.\d+)\s*\(\d+\)\s*(\d\.\d+)\s*\(\d+\)\s*(\d\.\d+)\s\(\d+\)\s*" % key
        match_obj = re.search( regex, output )
        if match_obj:
            max_time = float( match_obj.group( 3 ) )
        else:
            sys.exit( "Could not find regex in string." )

        return max_time
# ==============================================================================
def _parse_options():
    '''Parse input options.'''
    import optparse, sys

    usage = "usage: %prog filename"

    parser = optparse.OptionParser( usage = usage )

    (options, args) = parser.parse_args()

    #if not args  or  len(args) != 1:
        #parser.print_help()
        #sys.exit( "\nProvide a file to be split." )

    return args[0]
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
        sys.exit( "\nERROR: The command \n\n%s\n\nreturned a nonzero " \
                  "exit status. The error message is \n\n%s\n\n" \
                  "Abort.\n" % \
                  ( command, process.stderr.read()[:-1] )
                )
    return output
# ==============================================================================
if __name__ == '__main__':
    STATUS = _main()
    sys.exit( STATUS )
# ==============================================================================
