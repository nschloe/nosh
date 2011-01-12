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

    # speedup
    speedup = min_vals[0] / min_vals * num_procs[0]
    pp.plot( num_procs, speedup, 'ok' )
    pp.plot( [0,50], [0,50], '-k' )
    pp.title( "Speedup" )
    pp.xlabel( "Number of processors" )
    pp.ylabel( "speedup" )
    pp.show()
    #matplotlib2tikz.save( "speedup.tikz" )

    # efficiency
    efficiency = speedup / num_procs
    pp.plot( num_procs,  efficiency, 'ok' )
    pp.plot( [0,50], [1,1], '-k' )
    pp.title( "Efficiency" )
    pp.show()
    #matplotlib2tikz.save( "efficiency.tikz" )

    # serial fraction
    serial_fraction = (1.0/speedup - 1.0/num_procs) / (1.0 - 1.0/num_procs)
    pp.plot( num_procs,  serial_fraction, 'ok' )
    #pp.plot( [0,50], [1,1], '-k' )
    pp.title( "Serial fraction" )
    pp.show()

    return
# ==============================================================================
def _read_data( filename ):

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

    num_entries = len( num_procs )

    # convert them into int
    num_procs = np.array( num_procs, dtype = int )

    max_numrows = 1000
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
