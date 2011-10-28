#!/usr/bin/env python
'''Plots a file with scalability timings by taking the minimal values
across a column for all columns of the data file, and plot the ratio
of the value and the value of the first column.
Strong scaling.
'''
# ==============================================================================
import numpy as np
import math
import sys
# ==============================================================================
def _main():
    import matplotlib.pyplot as pp
#    import matplotlib2tikz

    filenames, is_tikz = _parse_options()

    num_procs = []
    labels    = []
    min_vals  = []
    max_procs = []
    for filename in filenames:
        label, num_p, mv = _read_data( filename )
        labels.append( label )
        num_procs.append( num_p )
        min_vals.append( mv )
        max_procs.append( max( num_p ) )

    # get the overall maximum
    max_procs = max( max_procs )

    # plot the data
    marker_styles = [ '+', '*', '1', '.', ',', '2', '3', '4', '<', '>', 'D',
                      'H', '^', '_', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|'
                    ]
    # --------------------------------------------------------------------------
    # Plot the times.
    # First plot a reference line of ideal speedup.
    # Center this line at the log-average of all other plots.
    log10_average = 0
    for min_val in min_vals:
        log10_average = log10_average + math.log10( min_val[0] )
    log10_average = log10_average / len(min_vals)
    proc_range = range( 1, max_procs+1 )
    pp.plot( proc_range,
             10**log10_average/np.array(proc_range),
             '-k',
             linewidth = 2,
             label = "ideal speedup"
           )
    pp.title( "Epetra timing for shared-memory, MPI" )
    pp.xlabel( "Number of processes" )
    pp.ylabel( "time (s)" )
    pp.xlim( 1, max_procs+1 )
    k = 0
    for min_val, num_proc, label  in zip( min_vals, num_procs, labels ):
        pp.loglog( num_proc, min_val,
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
    pp.xlabel( "Number of processes" )
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
    pp.xlabel( "Number of processes" )
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

    # read the
    header = data_reader.next()[0]
    # strip the initial "#" and split by ';'
    header = header.lstrip( '#' ).split( ';' )
    label = header[0].strip()

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

    return label, num_procs, min_values
# ==============================================================================
def _parse_options():
    '''Parse input options.'''
    import argparse

    parser = argparse.ArgumentParser( description = 'Plot scalability data as '
                                            + 'produced by scalability_timer.' )

    parser.add_argument( 'filenames',
                         metavar = 'FILE',
                         nargs   = '+',
                         type    = str,
                         help    = 'file to be read from'
                       )

    parser.add_argument( '--tikz',
                         dest    = 'is_tikz',
                         action  = 'store_const',
                         const   = True,
                         default = False,
                         help='produce TikZ/Pgfplots output (default: show on screen)')

    args = parser.parse_args()

    return args.filenames, args.is_tikz
# ==============================================================================
if __name__ == "__main__":
    _main()
# ==============================================================================
