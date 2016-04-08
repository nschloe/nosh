# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
from .helpers import extract_c_expression, templates_dir


def get_code_expression(name, expr):
    x = sympy.MatrixSymbol('x', 3, 1)

    result = self.expr.eval(x)

    expr_arguments = set([sympy.Symbol('x')])
    try:
        free_vars = result.free_variables
    except AttributeError:
        free_vars = set([])
    unused_arguments = expr_arguments - free_vars

    # template substitution
    with open(os.path.join(templates_dir, 'expression.tpl'), 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': self.name.lower(),  # class names are lowercase
            'degree': self.expr.degree,
            'eval': extract_c_expression(result),
            'unused_args': '\n'.join(
                ('(void) %s;' % name) for name in unused_arguments
                ),
            })

    return code, set()
