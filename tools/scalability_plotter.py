#!/usr/bin/env python
'''Plots a file with scalability timings by taking the minimal values
across a column for all columns of the data file, and plot the ratio
of the value and the value of the first column.
Strong scaling.
'''
# ==============================================================================
import numpy as np
import sys
# ==============================================================================
def _main():
    import matplotlib.pyplot as pp
    import matplotlib2tikz

    filenames, is_tikz = _parse_options()

    num_procs = []
    min_vals = []
    max_procs = []
    for filename in filenames:
        num_p, mv = _read_data( filename )
        num_procs.append( num_p )
        min_vals.append( mv )
        max_procs.append( max( num_p ) )

    # get the overall maximum
    max_procs = max( max_procs )

    labels = [ "Vector::MeanValue",
             "Vector::MinValue",
             "Vector::MaxValue",
             "Vector::Norm1",
             "Vector::Norm2",
             "Vector::NormInf",
             "Vector::Scale",
             "Vector::Dot",
             "Vector::Multiply",
             "Vector::Update",
             "CrsMatrix::Norm1",
             "CrsMatrix::NormInf",
             "CrsMatrix::NormFrobenius",
             "CrsMatrix::Scale",
             "CrsMatrix::LeftScale",
             "CrsMatrix::RightScale",
             "CrsMatrix::Apply" ]

    # plot the data
    marker_styles = [ '+', '*', '1', '.', ',', '2', '3', '4', '<', '>', 'D',
                      'H', '^', '_', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|'
                    ]
    # --------------------------------------------------------------------------
    # times
    proc_range = range( 1, max_procs+1 )
    pp.plot( proc_range, 1.0/np.array(proc_range), '-k', label="ideal speedup" )
    pp.title( "Epetra timing for shared-memory, MPI" )
    pp.xlabel( "Number of processes" )
    pp.xlim( 0, max_procs+1 )
    k = 0
    for min_val, num_proc, label  in zip( min_vals, num_procs, labels ):
        pp.semilogy( num_proc, min_val,
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
def _read_data( filename ):
    '''Reads the CSV data from a file and returns the processed data.'''
    import csv

    # initialize csv reader
    data_reader = csv.reader( open(filename, 'r'),
                              delimiter='\t',
                              quotechar='|',
                              skipinitialspace = True
                            )

    # skip the header
    data_reader.next()

    # Read num_procs and remove empty entries by list comprehension
    num_procs = [x for x in data_reader.next() if x]

    # convert them into int
    num_procs = np.array( num_procs, dtype = int )

    # Allocate space for the values.
    num_entries = len( num_procs )
    max_numrows = 10000
    values = np.empty( (max_numrows, num_entries), dtype = float )

    # read all the data into an array
    try:
        k = 0
        for row in data_reader:
            data = [x for x in row if x] # filter out empty values
            values[ k, :len(data) ] = data
            if len( data ) < num_entries:
                # pad with inf
                values[ k, len(data): ] = np.inf
            k += 1
    except csv.Error, e:
        sys.exit( 'file %s, line %d: %s' % (filename, data_reader.line_num, e) )

    # trim the values
    values = values[ :k, : ]

    # take the minimum along all rows
    min_values = np.amin( values, 0 )

    return num_procs, min_values
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
