#!/usr/bin/env python
'''Performs scalability runs, i.e., runs a given code through mpiexec -n X
with various X, parses the output and stores it in a file. Typically, one
would store timings then of course.'''
# ==============================================================================
import sys, subprocess, re, datetime, os
# ==============================================================================
def _main():

    # define the command to execute
    cmd_serial = "/home/nschloe/ginla/build/mpi/examples/linear-solve-test/keo-belos.exe --input=./cube.e"
    cmd_parallel = "mpiexec -n NUMPROCS /home/nschloe/ginla/build/mpi/examples/linear-solve-test/keo-belos.exe --input=./cube-balanced.par"
    # specify the keys to look for in timing outputs
    keys = [ "Belos: PseudoBlockCGSolMgr total solve time" ]
    # specify the number of procs to run this problem on
    num_procs = range( 1, 49 )
    # specify how many runs to do for each individual num_proc
    num_runs = 10
    # provide a comment to be included in the logs
    comment = "Timings of Belos solves with a Jacobian system, preconditioned" \
            + " with KEO (solved with ML)"

    # --------------------------------------------------------------------------
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

    # write header to timing files
    hostname = os.uname()[1]
    time = datetime.datetime.now().isoformat(' ')
    for handle, key in zip( timingfile_handles, keys ):
        header = "# %s; %s; %s\n" % ( key, hostname, time )
        handle.write( header )

    # write header to stdout file
    header = "# %s, %s\n" % ( hostname, time )
    outputfile_handle.write( header )

    # write number of processors to timing files
    for handle in timingfile_handles:
        for np in num_procs:
            handle.write( "%d\t\t" % np )
        handle.write( "\n" )

    # perform the test runs
    for k in xrange( 0, num_runs ):
        # ----------------------------------------------------------------------
        for np in num_procs:
            # run the test
            if np == 1:
                output, max_times = _testrun( cmd_serial, np, keys )
            elif np > 1:
                cmd = cmd_parallel.replace( "NUMPROCS", str(np) )
                output, max_times = _testrun( cmd, np, keys )
            else:
                sys.exit( "Illegal number of processors \"%d\"." % np )

            # write timing data
            for handle, max_time in zip( timingfile_handles, max_times ):
                handle.write( "%e\t" % max_time )

            # write stdout
            outputfile_handle.write( 2*(80*"#" + "\n")
                                     + "Run on %d cores, loop no. %d\n" \
                                       % ( np, k )
                                     + 80*"=" + "\n\n" )
            outputfile_handle.write( output )
            outputfile_handle.write( "\n" )
        # ----------------------------------------------------------------------
        for handle in timingfile_handles:
            handle.write( "\n" )

    outputfile_handle.write( 2*(80*"#" + "\n")
                             + "EOF" )

    # close all files
    outputfile_handle.close()
    for handle in timingfile_handles:
        handle.close()

    return
# ==============================================================================
def _testrun( cmd, num_procs, keys ):

    # regular expression for a floating point number
    fp_regex = '\d+\.?\d*(?:[eE][+-]\d+)?'

    regexs = []
    if num_procs == 1:
        # "Belos: PseudoBlockCGSolMgr total solve time    0.3433 (1)  "
        for key in keys:
            regexs.append( "%s\s*(%s)" % ( key, fp_regex ) )
    elif num_procs > 1:
        # "Some random keyword    0.3433 (1)  0.3434 (1)  0.03435 (1)   "
        # Store the *last* timing (i.e., the max across all processors)
        for key in keys:
            regexs.append( "%s\s*%s\s*\(\d+\)\s*%s\s*\(\d+\)\s*(%s)\s*\(\d+\)" \
                           % (key, fp_regex, fp_regex, fp_regex)
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
