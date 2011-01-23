#!/usr/bin/env python
'''Plots a file with scalability timings by taking the minimal values
across a column for all columns of the data file, and plot the ratio
of the value and the value of the first column.
Strong scaling.
'''
# ==============================================================================
import numpy as np
# ==============================================================================
def _main():
    import matplotlib.pyplot as pp
    import matplotlib2tikz

    filename, is_tikz = _parse_options()

    num_procs, min_vals = _read_data( filename )

    max_procs = max( num_procs )

    # plot the data

    # speedup
    speedup = min_vals[0] / min_vals * num_procs[0]
    pp.plot( num_procs, speedup, 'ok', label="speedup" )
    pp.plot( [0,max_procs+1], [0,max_procs+1], '-k', label="ideal" )
    pp.title( "Speedup" )
    pp.xlabel( "Number of processors" )
    pp.ylabel( "speedup" )
    pp.xlim( 0, max_procs+1 )
    pp.ylim( 0, max_procs+1 )
    pp.legend()
    if is_tikz:
        matplotlib2tikz.save( "speedup.tikz" )
    else:
        pp.show()

    # efficiency
    efficiency = speedup / num_procs
    pp.plot( num_procs,  efficiency, 'ok', label="efficiency" )
    pp.plot( [0,max_procs+1], [1,1], '-k', label="ideal" )
    pp.title( "Efficiency" )
    pp.xlim( 0, max_procs+1 )
    pp.ylim( 0, 1.1 )
    pp.legend()
    if is_tikz:
        matplotlib2tikz.save( "efficiency.tikz" )
    else:
        pp.show()

    # serial fraction
    serial_fraction = (1.0/speedup - 1.0/num_procs) / (1.0 - 1.0/num_procs)
    pp.plot( num_procs,  serial_fraction, 'ok' )
    pp.title( "Serial fraction" )
    if is_tikz:
        matplotlib2tikz.save( "serial-fraction.tikz" )
    else:
        pp.show()

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
            values[k,:] = filter( None, row )
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

    usage = "usage: %prog filename"

    parser = optparse.OptionParser( usage = usage )

    parser.add_option( "-t", "--tikz",
                       action  = "store_true",
                       dest    = "is_tikz",
                       default = False,
                       help    = "Plot the figures as TikZ files",
                       metavar = "TIKZFILE"
                     )

    (options, args) = parser.parse_args()

    if not args  or  len(args) != 1:
        parser.print_help()
        sys.exit( "\nProvide a file to be read." )

    return args[0], options.is_tikz
# ==============================================================================
if __name__ == "__main__":
    _main()
# ==============================================================================
