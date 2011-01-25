#!/usr/bin/env python
'''Plots a scalability timing data that is read from Trilinos output files
with a very specific format. Intended to be used with output from CalcUA.
'''
# ==============================================================================
import numpy as np
import sys
# ==============================================================================
def _main():
    import matplotlib.pyplot as pp
    import matplotlib2tikz

    filenames, is_tikz = _parse_options()

    # determine labels and matching regex keys
    labels = [ 'FECrsMatrix::Apply()',
               'Vector::Norm2()',
               'Vector::Dot()',
               'Vector::Multiply()',
             ]
    keys = [ "Matrix-vector multiplication",
             "2-norm calculation",
             "inner product",
             "element-wise multiplication"
           ]

    # get the data
    data = _parse_data( filenames, labels, keys )

    # reduce times to minimal value
    for key, value in data.iteritems():
        for numprocs in value['times'].keys():
            value['times'][numprocs] = min( value['times'][numprocs] )

    # plot the data
    marker_styles = [ '+', '*', '1', '.', ',', '2', '3', '4', '<', '>', 'D',
                      'H', '^', '_', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|'
                    ]
    # --------------------------------------------------------------------------
    # times
    #proc_range = range( 1,  )
    #pp.plot( proc_range, 1.0/np.array(proc_range), '-k', label="ideal speedup" )
    pp.title( "Epetra timing for shared-memory, MPI" )
    pp.xlabel( "Number of processes" )
    #pp.xlim( 0, max_procs+1 )
    k = 0

    for label, value in data.iteritems():
        # sort by num_procs
        numprocs = value['times'].keys()
        numprocs.sort()
        minvals = map( value['times'].get, numprocs )

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
def _parse_data( filenames, labels, keys ):
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
'''
    import re

    # initialize the data structure
    data = {}
    for (label,key) in zip(labels,keys):
        data[ label ] = {}
        data[ label ][ 'key' ]   = key
        data[ label ][ 'times' ] = {}

    # loop over the files
    for filename in filenames:
        # open the file
        file_handle = open( filename, 'r' )

        # prepare for iterating through the file
        content = file_handle.read()

        # generate regex for number of processes
        regex0 = "# (\d+) processes"

        # parse for num_procs

        match_obj = re.search( regex0, content )
        if match_obj:
            num_procs = int( match_obj.group( 1 ) )
        else:
            error_message = "Could not find the regex \n\n%s\n\n " \
                            "in the string \n\n%s\n\n." \
                          % ( regex, repr(output) )
            sys.exit( error_message )

        # generate regular expressions for keys
        regexs = []
        if num_procs == 1:
            # "Belos: PseudoBlockCGSolMgr total solve time    0.3433 (1)  "
            for key in keys:
                regexs.append( "%s\s*(\d+\.?\d*)" % key )
        elif num_procs > 1:
            # "Some random keyword    0.3433 (1)  0.3434 (1)  0.03435 (1)   "
            # Store the *last* timing (i.e., the max across all processors)

            # regular expression for a floating point number:
            fp_regex = '[-+]?\d*\.?\d*(?:e[-+]?\d+)?'
            #fp_regex = '\d*\.?\d*'
            for key in keys:
                regexs.append( "%s\s*%s\s*\(\d+\)\s*%s\s*\(\d+\)\s*(%s)\s*\(\d+\)" \
                              % (key, fp_regex, fp_regex, fp_regex)
                            )
        else:
            sys.exit( "Illegal number of processors \"%d\"." % num_procs )

        # parse the output for the keys
        for regex, label in zip(regexs, labels):
            # initialize list
            if not num_procs in data[label]['times'].keys():
                data[label]['times'][num_procs] = []
            match_obj = re.search( regex, content )
            if match_obj:
                # append time value
                time = float( match_obj.group( 1 ) )
                data[label]['times'][num_procs].append( time )
            else:
                error_message = "Could not find the regex \n\n%s\n\n " \
                                "in the string \n\n%s\n\n." \
                              % ( regex, repr(content) )
                sys.exit( error_message )

    return data
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
