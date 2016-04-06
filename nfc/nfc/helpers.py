# -*- coding: utf-8 -*-
#
import os
import re
import subprocess
import sys

templates_dir = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '..',
    'templates'
    )


def extract_c_expression(expr):
    from sympy.utilities.codegen import codegen
    [(c_name, c_code), (h_name, c_header)] = codegen(("f", expr), "C")
    res = re.search("f_result = (.*);", c_code)
    return res.group(1)


def run(command):
    """Runs a given command on the command line and returns its output.
    """
    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        close_fds=True
        )
    output = process.stdout.read()[:-1]
    ret = process.wait()

    if ret != 0:
        sys.exit(
            "\nERROR: The command \n\n%s\n\nreturned a nonzero "
            "exit status. The error message is \n\n%s\n\n"
            "Abort.\n"
            % (command, process.stderr.read()[:-1])
            )
    return output
