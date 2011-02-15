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

    # Get the data.
    # The data is a nested dict of dicts of dicts, in the order
    #    data[ host ][ function ] = [ [num_procs],
    #                                 [[time_data]]
    #                               ]
    data, hosts = _parse_data( filenames )

    if len( data ) == 0:
        sys.exit( '\n\tNo suitable data found at all. Abort.\n' )

    # add minimal value
    for host, data0 in data.iteritems():
        for function_name, data1  in data0.iteritems():
            for num_procs, data2 in data1.iteritems():
                data[host][function_name][num_procs] = [
                                          data[host][function_name][num_procs],
                                     min( data[host][function_name][num_procs] )
                                                       ]

    # Now the data is of the form
    #
    # The data is a nested dict of dicts of dicts, in the order
    #    data[ host ][ function ][ num_procs ] = [ [time_data],
    #                                              min( time_data )
    #                                            ]

    # Sort the entries by number of processes.
    # This is to make sure that when the data is plotted further below, and
    # lines are used, only subsequent numbers of processors are connected.
    for host, data0 in data.iteritems():
        for function_name, data1  in data0.iteritems():
            # Sort by num_procs.
            # See <http://code.activestate.com/recipes/52306/>.
            numprocs_list = data1.keys()
            numprocs_list.sort()
            print numprocs_list
            print map( data1.get, numprocs_list )
            data[host][function_name] = ( numprocs_list,
                                          map( data1.get, numprocs_list )
                                        )

    # get smallest and largest number of cores used
    min_procs = np.inf
    max_procs = -np.inf
    for host, data0 in data.iteritems():
        for function_name, data1  in data0.iteritems():
            min_procs = min( min_procs, min(data1[0]) )
            max_procs = max( max_procs, max(data1[0]) )

    # set marker styles
    marker_styles = [ '+', '*', '1', '.', ',', '2', '3', '4', '<', '>', 'D',
                      'H', '^', '_', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|'
                    ]
    # --------------------------------------------------------------------------
    # Plot the raw times.
    pp.title( ', '.join( hosts ) ) # use hostnames as title
    pp.xlabel( "Number of processes" )
    pp.xlim( 0, max_procs+1 )

    # plot ideal speedup
    proc_range = range( 1, max_procs+1 )
    pp.plot( proc_range,
             1.0 / np.array(proc_range),
             '-k',
             linewidth = 2,
             label = "ideal speedup"
           )

    k = 0
    for host, data0 in data.iteritems():
        for function_name, data1  in data0.iteritems():
            pp.semilogy( data1[0], data1[1],
                         linestyle = '-',
                         marker    = marker_styles[ k%len(marker_styles) ],
                         label     = function_name
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

    for host, data0 in data.iteritems():
        for function_name, data1  in data0.iteritems():
            for num_procs, data2 in data1.iteritems():
                speedup = minvals[0] / minvals * numprocs[0]
                data[host][function_name][num_procs].append (

                )


    for host, data0 in data.iteritems():
        for function_name, data1 in data0.iteritems():
            data[host][function_name]
            pp.semilogy( data1[0], data1[1],
                         linestyle = '-',
                         marker    = marker_styles[ k%len(marker_styles) ],
                         label     = function_name
                       )
            k += 1

    for label, value in data.iteritems():
        # sort by num_procs
        numprocs = value.keys()
        numprocs.sort()
        minvals = np.array( map( value.get, numprocs ) )

        speedup = minvals[0] / minvals * numprocs[0]
        speedups.append( speedup )
        pp.plot( numprocs, speedup,
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
    for label, value in data.iteritems():
        # sort by num_procs
        numprocs = value.keys()
        numprocs.sort()
        minvals = np.array( map( value.get, numprocs ) )

    for speedup, num_proc, label  in zip( speedups, num_procs, labels ):
        efficiency = speedup / num_proc
        pp.plot( numprocs, efficiency,
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

# arch: <some architecture specs>
# numprocs: 8
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

    # initialize the data structure
    data = {}
    hosts = []

    # loop over the files
    for filename in filenames:
        # open the file
        file_handle = open( filename, 'r' )

        # prepare for iterating through the file
        content = file_handle.read()

        # parse for system architecture
        arch = get_arch( content )
        if not arch:
            print "System architecture spec not found: Skipping file %s." \
                  % filename
            continue

        # parse for num_procs
        num_procs = get_numprocs( content )
        if not num_procs:
            print "Numprocs not found: Skipping file %s." % filename
            continue

        # generate regular expressions for keys
        regex = generate_regex( num_procs )

        try:
            # Make ^ and $ match the beginning and end of a line, respectively.
            match_obj_iter = re.finditer( regex, content, re.MULTILINE )
        except:
            print "Invalid file \"%s\". Skip." % filename
            continue

        # parse for hosts
        host = get_host( content )
        if not host:
            print "Host name not found: Skipping file %s." \
                  % filename
            continue
        # add to host list
        if host not in hosts:
            hosts.append( host )

        for match_obj in match_obj_iter:
            if match_obj:
                # get key and time value
                func_name = match_obj.group( 1 ).strip()
                time = float( match_obj.group( 2 ) )

                # add data in the dictionary of dictionaries
                if host not in data.keys():
                    data[ host ] = {}
                if func_name not in data[host].keys():
                    data[ host ][ func_name ] = [ [], [] ]

                # find num_procs in the list
                num_procs_list = data[host][func_name][ 0 ]
                num_procs_pos = [i for i,x in enumerate(num_procs_list) if x == num_procs]

                if len( num_procs_pos ) == 0:
                    # insert num_procs

                if num_procs not in data[host][func_name][ 0 ]:
                    data[ host ][ func_name ][ 0 ].append( num_procs )

                # Finally add the time value at the index where num_procs is
                # located
                data[ host ][ func_name ][ num_procs ][ 1 ][ ].append( time )

            else:
                error_message = "Could not find the regex \n\n%s\n\n " \
                                "in the string \n\n%s\n\n." \
                              % ( regex, repr(content) )
                sys.exit( error_message )

    return data, hosts
# ==============================================================================
def get_numprocs( content ):
    '''Parse a string for "# numprocs: 8" and extract the number of procs.'''
    regex = "^# numprocs: (\d+)$"

    match_obj = re.search( regex, content, re.MULTILINE )
    if match_obj:
        return int( match_obj.group( 1 ) )
    else:
        return None
# ==============================================================================
def get_arch( content ):
    '''Parse a string for "# 8 processes" and extract the number of procs.'''
    regex = "^# arch: (.*)$"

    match_obj = re.search( regex, content, re.MULTILINE )
    if match_obj:
        return match_obj.group( 1 )
    else:
        return None
# ==============================================================================
def get_host( content ):
    '''Parse a string for "# host: turing.'''
    regex = "^# host: (.*)$"

    match_obj = re.search( regex, content, re.MULTILINE )
    if match_obj:
        return match_obj.group( 1 )
    else:
        return None
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
