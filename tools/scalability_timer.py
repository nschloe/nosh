#!/usr/bin/env python
'''Performs scalability runs, i.e., runs a given code through mpiexec -n X
with various X, parses the output and stores it in a file. Typically, one
would store timings then of course.'''
# ==============================================================================
import sys, subprocess, re, datetime, os
# ==============================================================================
def _main():

    keys = [ "Matrix-vector multiplication",
             "2-norm calculation",
             "inner product",
             "element-wise multiplication"
           ]
    comments = [ "Timings of one matrix-vector product with KEO, Cube 3D",
                 "Timings of one 2-norm calculation, Cube 3D",
                 "Timings of one inner product calculation, Cube 3D",
                 "Timings of one element-wise vector-multiplication, Cube 3D"
               ]

    # set the data file
    timingfile_handles = []
    for k in xrange( len(keys) ):
        timing_file = "scaling-test-%d.dat" % k
        bufsize = 0 # write out the data immediately
        timingfile_handles.append( open( timing_file, "w", bufsize ) )

    # set the file for stdout
    output_file = "output.log"
    bufsize = 1
    outputfile_handle = open( output_file, "w", bufsize )

    # run over the number of procs
    min_numprocs = 1
    max_numprocs = 48

    num_runs = 10000

    #comment = "Timings of Belos solves with a Jacobian system, preconditioned
               #with KEO (solved with ML)"
    #comment = "Timings of one full LOCA step"

    # write header to timing files
    hostname = os.uname()[1]
    time = datetime.datetime.now().isoformat(' ')
    k = 0
    for handle in timingfile_handles:
        header = "# %s, %s, %s\n" % ( comments[k], hostname, time )
        handle.write( header )
        k += 1

    # write header to stdout file
    header = "# %s, %s\n" % ( hostname, time )
    outputfile_handle.write( header )

    # write number of processors to timingfiles
    for handle in timingfile_handles:
        for num_procs in xrange( min_numprocs, max_numprocs+1 ):
            handle.write( "%d\t\t" % num_procs )
        handle.write( "\n" )

    # perform the test runs
    for k in xrange( 0, num_runs ):
        # ----------------------------------------------------------------------
        for num_procs in xrange( min_numprocs, max_numprocs+1 ):
            # run the test
            output, max_times = _testrun( num_procs, keys )

            # write timing data
            for kk in xrange( len(max_times) ):
                timingfile_handles[kk].write( "%e\t" % max_times[kk] )

            # write stdout
            outputfile_handle.write( 2*(80*"#" + "\n")
                                     + "Run on %d cores, loop no. %d\n" \
                                       % ( num_procs, k )
                                     + 80*"=" + "\n\n" )
            outputfile_handle.write( output )
            outputfile_handle.write( "\n" )
        # ----------------------------------------------------------------------
        for handle in timingfile_handles:
            handle.write( "\n" )

    # close all files
    outputfile_handle.close()
    for handle in timingfile_handles:
        handle.close()

    return
# ==============================================================================
def _testrun( num_procs, keys ):

    test_exe = "/home/nschloe/ginla/build/mpi/packages/ginla-fvm/examples/linear-solve-test/keo-matvec.exe"
    basename = "cube"

    #test_exe = "/home/nschloe/ginla/build/mpi/packages/ginla-fvm/examples/linear-solve-test/keo-belos.exe"
    #basename = "cutcircle1000"
    #key = "Belos: PseudoBlockCGSolMgr total solve time"

    #test_exe = "/home/nschloe/ginla/build/mpi/packages/ginla-fvm/examples/loca-driver/loca-driver.exe"
    #basename = "cutcircle1000"
    #key = "Belos: PseudoBlockCGSolMgr total solve time"

    #test_exe = "/home/nico/ginla/build/mpi/packages/ginla-fvm/examples/linear-solve-test/keo-belos.exe"
    #basename = "cutcircle300"
    #key = "Belos: PseudoBlockCGSolMgr total solve time"

    if num_procs == 1:
        options = "--input=%s.e" % basename
    elif num_procs > 1:
        options = "--input=%s-balanced.par" % basename
    else:
        sys.exit( "Illegal number of processors \"%d\"." % num_procs )

    #test_exe = "/home/nschloe/ginla/build/mpi/packages/ginla-fvm/examples/loca-driver/loca-driver-fvm.exe"
    #options = "--xml-input-file=./conf.xml"
    #key = "LOCA runtime"

    if num_procs == 1:
        cmd = "%s %s" % ( test_exe, options )
        # "Belos: PseudoBlockCGSolMgr total solve time    0.3433 (1)  "
        regexs = []
        for key in keys:
            regexs.append( "%s\s*(\d+\.?\d*)" % key )
    elif num_procs > 1:
        cmd = "mpiexec -n %d %s %s" % ( num_procs, test_exe, options )
        # "Some random keyword    0.3433 (1)  0.3434 (1)  0.03435 (1)   "
        # Store the *last* timing (i.e., the max across all processors)
        regexs = []
        decimal_regex = '\d+\.?\d*' # regular expression for a decimal number
        for key in keys:
            regexs.append( "%s\s*%s\s*\(\d+\)\s*%s\s*\(\d+\)\s*(%s)\s*\(\d+\)" \
                           % (key, decimal_regex, decimal_regex, decimal_regex)
                         )
    else:
        sys.exit( "Illegal number of processors \"%d\"." % num_procs )

    # run the test
    output = _run( cmd )

    max_times = []
    for regex in regexs:
        # Parse the output, get the maximum execution time over all processes.
        # Find something of the kind
        match_obj = re.search( regex, output )
        if match_obj:
            max_times.append( float( match_obj.group( 1 ) ) )
        else:
            error_message = "Could not find the regex \n\n%s\n\n " \
                            "in the string \n\n%s\n\n." \
                          % ( regex, repr(output) )
            sys.exit( error_message )

    return output, max_times
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
