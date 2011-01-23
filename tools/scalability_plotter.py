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
        np, mv = _read_data( filename )
        num_procs.append( np )
        min_vals.append( mv )
        max_procs.append( max( np ) )

    # get the overall maximum
    max_procs = max( max_procs )

    labels = [ 'Apply()', 'Norm2()', 'Dot()', 'Multiply()' ]

    # plot the data
    marker_styles = [ '+', '*', '1', '.', ',', '2', '3', '4', '<', '>', 'D',
                      'H', '^', '_', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|'
                    ]
    # --------------------------------------------------------------------------
    # speedup
    pp.plot( [0,max_procs+1], [0,max_procs+1], '-k', label="ideal" )
    pp.title( "Speedup" )
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
    pp.plot( [0,max_procs+1], [1,1], '-k', label="ideal" )
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
    import csv

    # initialize csv reader
    data_reader = csv.reader( open(filename, 'r'),
                              delimiter='\t',
                              quotechar='|',
                              skipinitialspace = True
                            )

    # skip the header
    data_reader.next()

    # Read num_procs and remove empty entries by filter.
    num_procs = filter( None, data_reader.next() )
    # convert them into int
    num_procs = np.array( num_procs, dtype = int )

    # Allocate space for the values.
    num_entries = len( num_procs )
    max_numrows = 10000
    values = np.empty( (max_numrows,num_entries), dtype = float )

    # read all the data into an array
    try:
        k = 0
        for row in data_reader:
            # 'filter' out empty values
            data = filter( None, row )
            values[k,:len(data)] = data
            if len( data ) < num_entries:
                # pad with inf
                values[k,len(data):] = np.inf
            k += 1
    except csv.Error, e:
        sys.exit('file %s, line %d: %s' % (filename, reader.line_num, e))

    # trim the values
    values = values[:k,:]

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
