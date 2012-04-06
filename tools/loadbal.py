#!/usr/bin/env python
# ==============================================================================
import argparse, sys, os, subprocess, shutil
# ==============================================================================
def _main():
    args = _parse_options()

    # get basename of input file
    basename, extension = os.path.splitext( args.filename )

    # define output filename(s)
    output = "%s-balanced.nemI" % basename

    # slice it
    print "Cutting the input data into slices of ", args.proc_mesh, "."
    #for proc_mesh in args.proc_meshes:
    _slice(args.filename, output, args.proc_mesh, args.slice_method)

    return
# ==============================================================================
def _slice(filename, output, proc_mesh, slice_method):

    bin_dir = '/opt/trilinos/dev/master/gcc/4.6.1/bin/'
    #bin_dir = '/opt/trilinos/dev/master/openmpi/1.4.3/gcc/4.6.1/release/bin/'
    nemslice_command = bin_dir + "nem_slice"
    nemspread_command = bin_dir + "nem_spread"
    tmp_nemspreadinp = "nem_spread.inp"

    slice_command = "%s -v -o \"%s\" -e -m mesh=%dx%d -l %s \"%s\"" % \
                    ( nemslice_command, output, proc_mesh[0], proc_mesh[1], slice_method, filename )

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

    parser.add_argument('filename',
                        metavar='FILE',
                        type=str,
                        help='the file to read from')

    parser.add_argument('--process-mesh', '-p',
                        dest='proc_mesh',
                        metavar='PROCMESH',
                        type=int,
                        nargs=2,
                        required=True,
                        help='process mesh (e.g., -p 4 2)')

    parser.add_argument('--slice-method', '-m',
                        dest='slice_method',
                        type=str,
                        choices = ['inertial','spectral','multikl'],
                        default = 'multikl',
                        help='slicing method (default: multikl)')

    args = parser.parse_args()

    return args
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
