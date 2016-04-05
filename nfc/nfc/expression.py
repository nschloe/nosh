# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
from .helpers import extract_c_expression, templates_dir


class CodeExpression(object):
    def __init__(self, expr, name):
        self.expr = expr
        self.name = name

    def get_dependencies(self):
        return []

    def get_code(self):
        x = sympy.DeferredVector('x')

        result = self.expr.eval(x)

        # TODO
        # Check if any of the arguments is not used in the function.
        # (We'll declare them (void) to supress compiler warnings.)

        # template substitution
        with open(os.path.join(templates_dir, 'expression.tpl'), 'r') as f:
            src = Template(f.read())
            code = src.substitute({
                'name': self.name.lower(),  # class names are lowercase
                'degree': self.expr.degree,
                'eval': extract_c_expression(result)
                })

        return code
