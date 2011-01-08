#!/usr/bin/env python
# ==============================================================================
import sys, os, subprocess, shutil
# ==============================================================================
def _main():
        filename = _parse_options()

        # get basename of input file
        basename, extension = os.path.splitext( filename )

        # define output filename(s)
	output = "%s-balanced.nemI" % basename

        # slice it
        min_slices = 1
        max_slices = 48
        for num_slices in xrange( min_slices, max_slices+1 ):
            _slice( filename, output, num_slices )

	return
# ==============================================================================
def _slice( filename, output, num_slices ):

        nemslice_command = "nem_slice.exe"
        nemspread_command = "nem_spread.exe"
        slice_method = "spectral"

        slice_command = "%s -v -o \"%s\" -e -m mesh=1x%d -l %s \"%s\"" % \
                        ( nemslice_command, output, num_slices, slice_method, filename )

        _run( slice_command )

        # create a temporary directory in the current directory
        tmpdir = "tmp1"
        os.mkdir( tmpdir )

        # create a temporary file; nem_spread.inp
        file_handle = open( tmpdir + "/nem_spread.inp",
                            mode = "w" )
        #
        file_handle.write( "Input FEM file          = %s \
	LB file                 = %s \
	Debug                   = 1 \
	Restart Time list       = off \
	Reserve space           = nodal=1, elemental=0, global=0 \
	Parallel Disk Info = number=1 \
	Parallel file location = root=tmp,subdir=.. \
	" % ( filename, output ) )

        _run( nemspread_command )

        # remove the directories
	shutil.rmtree( tmpdir )

        return
# ==============================================================================
def _parse_options():
    '''Parse input options.'''
    import optparse, sys

    usage = "usage: %prog filename"

    parser = optparse.OptionParser( usage = usage )

    (options, args) = parser.parse_args()

    if not args  or  len(args) != 1:
        parser.print_help()
        sys.exit( "\nProvide a file to be split." )

    return args[0]
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
