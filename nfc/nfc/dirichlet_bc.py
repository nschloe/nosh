# -*- coding: utf-8 -*-
#
import nfl
import os
import sympy
from .helpers import templates_dir


class CodeDirichletBc(object):
    def __init__(self, dbc, name):
        assert(isinstance(dbc, nfl.DirichletBC))
        self.dbc = dbc
        self.name = name

    def get_dependencies(self):
        return []

    def get_code(self):
        x = sympy.DeferredVector('x')

        result0 = self.dbc.is_inside(x)
        try:
            is_x_used_inside = x in result0.free_symbols
        except AttributeError:
            is_x_used_inside = False

        result1 = self.dbc.eval(x)
        try:
            is_x_used_eval = x in result1.free_symbols
        except AttributeError:
            is_x_used_eval = False

        # template substitution
        with open(os.path.join(templates_dir, 'dirichlet_bc.tpl'), 'r') as f:
            src = Template(f.read())
            code = src.substitute({
                'name': self.name.lower(),  # class names are lowercase
                'inside_condition': extract_c_expression(result0),
                'inside_void': '' if is_x_used_inside else '(void) x;\n',
                'eval_return_value': extract_c_expression(result1),
                'eval_void': '' if is_x_used_eval else '(void) x;\n',
                })

        return code
