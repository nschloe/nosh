#!/usr/bin/env python
'''Plots a file with scalability timings by taking the minimal values
across a column for all columns of the data file, and plot the ratio
of the value and the value of the first column.
Strong scaling.
'''
# ==============================================================================
import csv
import numpy as np
import matplotlib.pyplot as pp
import matplotlib2tikz
# ==============================================================================
def _main():
    filename = _parse_options()

    num_procs, min_vals = _read_data( filename )

    # plot the data
    pp.plot( num_procs, min_vals[0] / min_vals * num_procs[0], 'ok' )
    pp.plot( [0,50], [0,50], '-k' )
    pp.title( "Speedup for matrix-vector product" )
    pp.xlabel( "Number of processors" )
    pp.ylabel( "speedup" )
    pp.show()
    #matplotlib2tikz.save( "speedup.tikz" )

    pp.plot( num_procs,  min_vals[0] / min_vals * num_procs[0] / num_procs, 'ok' )
    pp.plot( [0,50], [1,1], '-k' )
    pp.show()
    #matplotlib2tikz.save( "efficiency.tikz" )

    return
# ==============================================================================
def _read_data( filename ):

    # initialize csv reader
    data_reader = csv.reader( open(filename, 'r'),
                              delimiter='\t',
                              quotechar='|',
                              skipinitialspace = True
                            )

    # Read header and remove empty entries by filter.
    header = filter(None, data_reader.next())

    num_entries = len( header )

    # convert them into int
    num_procs = np.empty( num_entries, dtype = int )
    num_procs[:] = header

    max_numrows = 100
    values = np.empty( (max_numrows,num_entries), dtype = float )

    # read all the data into an array
    try:
        k = 0
        for row in data_reader:
            values[k,:] = row
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

    (options, args) = parser.parse_args()

    if not args  or  len(args) != 1:
        parser.print_help()
        sys.exit( "\nProvide a file to be read." )

    return args[0]
# ==============================================================================
if __name__ == "__main__":
    _main()
# ==============================================================================
