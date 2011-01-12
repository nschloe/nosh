#!/usr/bin/env python
'''Performs scalability runs, i.e., runs a given code through mpiexec -n X
with various X, parses the output and stores it in a file. Typically, one
would store timings then of course.'''
# ==============================================================================
import sys, subprocess, re
# ==============================================================================
def _main():
    filename = "scalingdata.dat"
    bufsize = 0 # write out the data immediately
    f = open( filename, "w", bufsize )

    # run over the number of procs
    min_numprocs = 1
    max_numprocs = 47

    num_runs = 10000

    # write header
    header_string = "#\t"
    for num_procs in xrange( min_numprocs, max_numprocs+1 ):
        header_string += "%d\t\t" % num_procs
    header_string += "\n"
    f.write( header_string )

    # perform the test runs
    for k in xrange( 0, num_runs ):
        for num_procs in xrange( min_numprocs, max_numprocs+1 ):
            max_time = _testrun( num_procs )
            f.write( "\t%e" % max_time )
        f.write( "\n" )

    f.close()

    return
# ==============================================================================
def _testrun( num_procs ):

    test_exe = "/home/nschloe/ginla/build/mpi/packages/ginla-fvm/examples/linear-solve-test/keo-belos.exe"
    basename = "cutcircle1000"

    key = "Belos: PseudoBlockCGSolMgr total solve time"

    if num_procs == 1:
        data_file = "%s.e" % basename
        cmd = "%s --input=%s" % ( test_exe, data_file )
        # "Belos: PseudoBlockCGSolMgr total solve time    0.3433 (1)  "
        regex = "%s\s*(\d+\.?\d*)" % key
    elif num_procs > 1:
        data_file = "%s-balanced.par" % basename
        cmd = "mpiexec -n %d %s --input=%s" % ( num_procs, test_exe, data_file )
        # "Belos: PseudoBlockCGSolMgr total solve time    0.3433 (1)  0.3434 (1)  0.03435 (1)   "
        # Enclose the *last* timing (i.e., the max across all processors) in parentheses
        regex = "%s\s*\d+\.?\d*\s*\(\d+\)\s*\d+\.?\d*\s*\(\d+\)\s*(\d+\.?\d*)\s*\(\d+\)" % key
    else:
        sys.exit( "Illegal number of processors \"%d\"." % num_procs )

    # run the test
    output = _run( cmd )

    # Parse the output, get the maximum execution time over all processes.
    # Find something of the kind
    match_obj = re.search( regex, output )
    if match_obj:
        max_time = float( match_obj.group( 1 ) )
    else:
        error_message = "Could not find the regex \n\n%s\n\n in the string \n\n%s\n\n." \
                      % ( regex, repr(output) )
        sys.exit( error_message )

    return max_time
# ==============================================================================
def _parse_options():
    '''Parse input options.'''
    import optparse

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
