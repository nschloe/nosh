#!/usr/bin/env python
'''Plots a scalability timing data that is read from Trilinos output files
with a very specific format. Intended to be used with output from CalcUA.
'''
# ==============================================================================
import numpy as np
import sys
import re
# ==============================================================================
def _main():
    import matplotlib.pyplot as pp
    import matplotlib2tikz

    filenames, is_tikz = _parse_options()

    # get the data
    data = _parse_data( filenames )

    # reduce times to minimal value
    for key, value in data.iteritems():
        for numprocs in value.keys():
            value[numprocs] = min( value[numprocs] )

    # plot the data
    marker_styles = [ '+', '*', '1', '.', ',', '2', '3', '4', '<', '>', 'D',
                      'H', '^', '_', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|'
                    ]
    # --------------------------------------------------------------------------
    # times
    #proc_range = range( 1,  )
    #pp.plot( proc_range, 1.0/np.array(proc_range), '-k', label="ideal speedup" )
    #pp.title( "Epetra timing for shared-memory, MPI" )
    pp.xlabel( "Number of processes" )
    #pp.xlim( 0, max_procs+1 )
    k = 0

    for label, value in data.iteritems():
        # sort by num_procs
        numprocs = value.keys()
        numprocs.sort()
        minvals = map( value.get, numprocs )

        pp.semilogy( numprocs, minvals,
                     linestyle = '-',
                     marker    = marker_styles[k],
                     label     = label
                   )
        k += 1
    pp.legend()
    if is_tikz:
        matplotlib2tikz.save( "timing.tikz" )
    else:
        pp.show()
    # --------------------------------------------------------------------------
    # speedup
    pp.plot( [0, max_procs+1], [0, max_procs+1], '-k', label="ideal" )
    pp.title( "Epetra speedup for shared-memory, MPI" )
    pp.xlabel( "Number of processes" )
    pp.xlim( 0, max_procs+1 )
    pp.ylim( 0, max_procs+1 )
    k = 0
    speedups = []
    for min_val, num_proc, label  in zip( min_vals, num_procs, labels ):
        speedup = min_val[0] / min_val * num_proc[0]
        speedups.append( speedup )
        pp.plot( num_proc, speedup,
                 linestyle = '-',
                 marker    = marker_styles[k],
                 label     = label
               )
        k += 1
    pp.legend( loc='upper left' )
    if is_tikz:
        matplotlib2tikz.save( "speedup.tikz" )
    else:
        pp.show()
    # --------------------------------------------------------------------------
    # efficiency
    pp.plot( [0, max_procs+1], [1, 1], '-k', label="ideal" )
    pp.title( "Efficiency" )
    pp.xlim( 0, max_procs+1 )
    pp.ylim( 0, 1.1 )
    k = 0
    for speedup, num_proc, label  in zip( speedups, num_procs, labels ):
        efficiency = speedup / num_proc
        pp.plot( num_proc, efficiency,
                 linestyle = '-',
                 marker    = marker_styles[k],
                 label     = label
               )
        k += 1
    pp.legend()
    if is_tikz:
        matplotlib2tikz.save( "efficiency.tikz" )
    else:
        pp.show()
    # --------------------------------------------------------------------------
    # serial fraction
    pp.title( "Serial fraction" )
    pp.xlim( 0, max_procs+1 )
    pp.ylim( 0, 1.1 )
    k = 0
    for speedup, num_proc, label  in zip( speedups, num_procs, labels ):
        serial_fraction = (1.0/speedup - 1.0/num_proc) / (1.0 - 1.0/num_proc)
        pp.plot( num_proc, serial_fraction,
                 linestyle = '-',
                 marker    = marker_styles[k],
                 label     = label
               )
        k += 1
    pp.legend()
    if is_tikz:
        matplotlib2tikz.save( "serial-fraction.tikz" )
    else:
        pp.show()
    # --------------------------------------------------------------------------

    return
# ==============================================================================
def _parse_data( filenames ):
    '''Parses Trilinos output of the form

# 8 processes
Restart Index set, reading solution time step: 1
Process 0 has 533732 nonzeros.
Process 6 has 538716 nonzeros.
Process 4 has 538952 nonzeros.
Process 7 has 547504 nonzeros.
Process 2 has 539656 nonzeros.
Process 3 has 539808 nonzeros.Process 1 has 537620 nonzeros.

Process 5 has 541652 nonzeros.
======================================================================================

                                 TimeMonitor Results

Timer Name                      Min over procs    Avg over procs    Max over procs
--------------------------------------------------------------------------------------
2-norm calculation              5.794e-05 (1)     0.0001656 (1)     0.00033 (1)
Data I/O                        4.871 (1)         5.27 (1)          5.573 (1)
FVM entities construction       0.1553 (1)        0.4671 (1)        0.8777 (1)
Graph construction              0.2614 (1)        0.2617 (1)        0.262 (1)
MVP construction                0.03962 (1)       0.05694 (1)       0.07495 (1)
Matrix construction             0.3576 (1)        0.3584 (1)        0.3603 (1)
Matrix-vector multiplication    0.002726 (1)      0.006423 (1)      0.00735 (1)
element-wise multiplication     0.0001121 (1)     0.0001621 (1)     0.000211 (1)
inner product                   0.005177 (1)      0.00518 (1)       0.005183 (1)
======================================================================================

The timer names are used as keys, and associated with the Max over procs,
split by the number of processes.
'''

    ## initialize the data structure
    data = {}
    #for (label,key) in zip(labels,keys):
        #data[ label ] = {}
        #data[ label ][ 'key' ]   = key
        #data[ label ][ 'times' ] = {}

    # loop over the files
    for filename in filenames:
        # open the file
        file_handle = open( filename, 'r' )

        # prepare for iterating through the file
        content = file_handle.read()

        # parse for num_procs
        num_procs = get_numprocs( content )

        # generate regular expressions for keys
        regex = generate_regex( num_procs )

        # Make ^ and $ match the beginning and end of a line, respectively.
        match_obj_iter = re.finditer( regex, content, re.MULTILINE )

        for match_obj in match_obj_iter:
            if match_obj:
                # get key and time value
                key = match_obj.group( 1 )
                time = float( match_obj.group( 2 ) )

                # add data structure
                if key not in data.keys():
                    data[key] = {}
                if num_procs not in data[key].keys():
                    data[key][num_procs] = []

                # append time value
                data[key][num_procs].append( time )
            else:
                error_message = "Could not find the regex \n\n%s\n\n " \
                                "in the string \n\n%s\n\n." \
                              % ( regex, repr(content) )
                sys.exit( error_message )

    return data
# ==============================================================================
def get_numprocs( content ):
    '''Parse a string for "# 8 processes" and extract the number of procs.'''
    regex = "# (\d+) processes"

    match_obj = re.search( regex, content )
    if match_obj:
        num_procs = int( match_obj.group( 1 ) )
    else:
        error_message = "Could not find the regex \n\n%s\n\n " \
                        "in the string \n\n%s\n\n." \
                      % ( regex, repr(content) )
        sys.exit( error_message )

    return num_procs
# ==============================================================================
def generate_regex( num_procs ):
    '''Generate regular expressions for the parsing of lines as
        Data I/O    4.871 (1)         5.27 (1)          5.573 (1)
    '''
    regexs = []

    # regular expression for a floating point number:
    fp_regex = '[-+]?\d*\.?\d*(?:e[-+]?\d+)?'

    # regex for " 0.3433 (1)"
    timing_regex       = "(?:%s)\s*\(\d+\)" % fp_regex
    timing_regex_store = "(%s)\s*\(\d+\)" % fp_regex

    if num_procs == 1:
        # "Belos: PseudoBlockCGSolMgr total solve time    0.3433 (1)  "
        regex = "^([^(?:%s)])*\s*(%s)\s*\(\d+\)$" \
                % (fp_regex, fp_regex)
    elif num_procs > 1:
        # "Some random keyword    0.3433 (1)  0.3434 (1)  0.03435 (1)   "
        # Store the *last* timing (i.e., the max across all processors)
        regex = "^(.*?)%s\s*%s\s*%s\s*$" \
                % ( timing_regex, timing_regex, timing_regex_store )
    else:
        sys.exit( "Illegal number of processors \"%d\"." % num_procs )

    return regex
# ==============================================================================
def _parse_options():
    '''Parse input options.'''
    import optparse

    usage = "usage: %prog datafile(s)"

    parser = optparse.OptionParser( usage = usage )

    parser.add_option( "-t", "--tikz",
                       action  = "store_true",
                       dest    = "is_tikz",
                       default = False,
                       help    = "Plot the figures as TikZ files",
                       metavar = "TIKZFILE"
                     )

    (options, args) = parser.parse_args()

    if not args:
        parser.print_help()
        sys.exit( "\nNo file(s) provided." )

    return args, options.is_tikz
# ==============================================================================
if __name__ == "__main__":
    _main()
# ==============================================================================
