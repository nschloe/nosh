#!/usr/bin/env python
# ==============================================================================
import re
# ==============================================================================
def _main():
    # Get the file handle from the command line, read the content into a string,
    # and close the file.
    loca_output_file = _parse_arguments()
    loca_file_content = loca_output_file.read()
    loca_output_file.close()

    # We're looking for the following string:
    #
    # The Belos solver of type "Belos::PseudoBlockCGSolMgr<...,double>{}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 36 iterations with total CPU time of 0.293983 sec
    # 
    # ************************************************************************
    # -- Nonlinear Solver Step 3 --
    # ||F|| = 3.651e-11  step = 1.000e+00  dx = 2.306e-06 (Converged!)
    # ************************************************************************
    real_number_regex    = '\d+\.\d+'
    # http://www.regular-expressions.info/floatingpoint.html
    #floating_point_regex = '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'
    floating_point_regex = '[0-9]+\.[0-9]+e[-+][0-9]+'
    # Create the monster of regular expression.
    regex = re.compile( r'''The Belos solver of type "Belos::PseudoBlockCGSolMgr<...,double>{}" returned a solve status of "SOLVE_STATUS_CONVERGED" in (\d+) iterations with total CPU time of %s sec

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*
-- Nonlinear Solver Step \d+ -- 
\|\|F\|\| = (%s)  step = (%s)  dx = (%s) \(Converged!\)
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*''' % ( real_number_regex, floating_point_regex, floating_point_regex, floating_point_regex ), re.MULTILINE )

    # Go through the contents of the file search for matches.
    print 0, 0 # For the 0th step (where the regex doesn't catch on).
    for ( k, match )  in enumerate( regex.finditer( loca_file_content ) ):
        print k+1, match.group( 1 )
    
    #print loca_output_file

    #for line in loca_output_file:
        #matches = regex.findall( line )
        #for match in matches:
            #print match

    #print re.search( regex, loca_output_file )

    return
# ==============================================================================
def _parse_arguments():
    '''Python 2.7 argument parser.'''

    import argparse

    parser = argparse.ArgumentParser( description='Extract information of the human-readable output of LOCA.' )

    parser.add_argument( 'loca_output_file',
                         type     = file,
                         help     = 'file containing the LOCA output'
                       )

    args = parser.parse_args()

    return args.loca_output_file
# ==============================================================================
if __name__ == "__main__":
    _main()
# ==============================================================================