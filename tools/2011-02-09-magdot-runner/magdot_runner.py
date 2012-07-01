#!/usr/bin/env python
'''Runs through different sizes of the retangular sample and the size of the
magnetic dot (that generates the potential), for each combination creating
a folder and start a continuation run in a detached session for it.
This is to
'''
# ==============================================================================
import sys, subprocess, os, shutil
from lxml import etree
import numpy as np
# ==============================================================================
def _main():

    exe = "/home/nschloe/ginla/build/mpi/packages/ginla-fvm/examples/loca-driver/loca-driver-fvm.exe"

    experiments_base_dir = '/home/nschloe/scratch/2011-02-14-magdot-enters-how'

    init_files = [ 'retangle-01.00.e',
                   'retangle-02.00.e',
                   'retangle-03.00.e',
                   'retangle-04.00.e',
                   'retangle-05.00.e',
                   'retangle-06.00.e',
                   'retangle-07.00.e',
                   'retangle-08.00.e',
                   'retangle-09.00.e',
                   'retangle-10.00.e'
                 ]

    radius_step = 0.5
    min_radius = radius_step
    max_radius = 5.0
    magDotRadii = np.arange( min_radius, max_radius+radius_step, radius_step )

    # read the XML config file
    conf_file_orig = 'conf.xml'
    conf_tree = etree.parse( conf_file_orig )
    # print etree.tostring(tree.getroot())

    for init_file in init_files:
        for radius in magDotRadii:
            basename = os.path.splitext(init_file)[0]

            # create a directory
            dir_name = os.path.join( experiments_base_dir,
                                     '%s-magdot-radius-%f' % ( basename, radius )
                                   )
            print "Creating directory '%s'..." % dir_name
            os.mkdir( dir_name )

            # Move the the init file there.
            shutil.copy( init_file, dir_name )

            # set the init file
            for element in conf_tree.iter( 'Parameter' ):
                if element.get( 'name' ) == 'State':
                    element.set( 'value', init_file )
                    break

            # set the magnetic field to 1.2345
            for element in conf_tree.iter( 'Parameter' ):
                if element.get( 'name' ) == 'magneticDotRadius':
                    element.set( 'value', '%f' % radius )
                    break

            # write the XML file
            conf_file = os.path.join( dir_name, 'conf.xml' )
            print conf_file
            f = open( conf_file, 'w' )
            f.write( etree.tostring( conf_tree.getroot() ) )
            f.close()

            # output file
            output_file = os.path.join( dir_name, 'output.log' )

            # start continuation as a background process
            opt = "--xml-input-file=%s" % conf_file
            # This following screen construct is to make sure we start the process
            # in a detached screen, and that the output is properly tee'ed to
            # a log file.
            # The "-d -m" options make sure that the session is started detached,
            # and the full command to be executed in the screen is given by an
            # sh string. "-c" for sh makes sure that the string is interpreted
            # as a command.
            cmd = 'screen -d -m sh -c "%s %s | tee %s"' % ( exe, opt, output_file )
            _run( cmd )


    return
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
