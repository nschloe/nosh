#!/usr/bin/env python
# ==============================================================================
import argparse, sys, os, subprocess, shutil
# ==============================================================================
def _main():
    filename = _parse_options()

    # get basename of input file
    basename, extension = os.path.splitext( filename )

    # define output filename(s)
    output = "%s-balanced.nemI" % basename

    # slice it
    slices_list = [12] #range( 2, 3 )
    #slices_list.extend( range(32, 33, 32) )
    print "Cutting the input data into slices of ", slices_list, "."
    for num_slices in slices_list:
        _slice( filename, output, num_slices )

    return
# ==============================================================================
def _slice( filename, output, num_slices ):

    bin_dir = "/opt/trilinos/dev/master/openmpi/1.4.3/gcc/4.6.1/release/bin/"
    nemslice_command = bin_dir + "nem_slice"
    nemspread_command = bin_dir + "nem_spread"
    tmp_nemspreadinp = "nem_spread.inp"
    slice_method = "inertial"
    #slice_method = "spectral"
    #slice_method = "multikl"

    slice_command = "%s -v -o \"%s\" -e -m mesh=1x%d -l %s \"%s\"" % \
                    ( nemslice_command, output, num_slices, slice_method, filename )

    _run( slice_command )

    # create a temporary directory in the current directory
    tmpdir = "tmp1"
    os.mkdir( tmpdir )

    # create a temporary file; nem_spread.inp
    file_handle = open( tmp_nemspreadinp,
                        mode = "w" )
    file_handle.write( "Input FEM file          = %s\n\
LB file                 = %s\n\
Debug                   = 1\n\
Restart Time list       = off\n\
Reserve space           = nodal=1, elemental=0, global=0\n\
Parallel Disk Info = number=1\n\
Parallel file location = root=tmp,subdir=.." % ( filename, output ) )
    file_handle.close()

    _run( nemspread_command )

    # clean up
    shutil.rmtree( tmpdir )
    os.remove( tmp_nemspreadinp )

    return
# ==============================================================================
def _parse_options():
    '''Parse input options.'''

    parser = argparse.ArgumentParser( description = 'Split an ExodusII file into loadbalanced chunks.' )

    parser.add_argument( 'filename',
                         metavar='FILE',
                         type=str,
                         nargs=1,
                         help='the file to read from')

    args = parser.parse_args()

    return args.filename[0]
# ==============================================================================
def _run( command ):
    """Runs a given command on the command line and returns its output.
    """
    print command

    process = subprocess.Popen( command,
                                shell     = True,
                                stdout    = subprocess.PIPE,
                                stderr    = subprocess.PIPE,
                                close_fds = True
                              )
    output = process.stdout.read()[:-1]
    ret = process.wait()

    if ret != 0:
        sys.exit( "\nERROR: The command \n\n%s\n\nreturned a nonzero " \
                  "exit status. The error message is \n\n%s\n\n" \
                  "Abort.\n" % \
                  ( command, process.stderr.read()[:-1] )
                )
    return output
# ==============================================================================
if __name__ == '__main__':
    STATUS = _main()
    sys.exit( STATUS )
# ==============================================================================
