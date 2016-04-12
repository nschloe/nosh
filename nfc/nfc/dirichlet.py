# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy

from .helpers import extract_c_expression, get_uuid, templates_dir
from .subdomain import *


class Dirichlet(object):
    def __init__(self, function, subdomains, is_matrix):
        self.class_name = 'dirichlet_' + get_uuid()
        self.dependencies = [SubdomainCode(sd) for sd in subdomains]
        self.function = function
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_obj(self):
        x = sympy.MatrixSymbol('x', 3, 1)
        result = self.function(x)
        try:
            is_x_used_eval = x in result.free_symbols
        except AttributeError:
            is_x_used_eval = False

        subdomain_ids = set([
            sd.__class__.__name__.lower() for sd in subdomains
            ])

        if len(subdomain_ids) == 0:
            # If nothing is specified, use the entire boundary
            subdomain_ids.add('boundary')

        init = '{%s}' % ', '.join(['"%s"' % s for s in subdomain_ids])

        # template substitution
        filename = os.path.join(templates_dir, 'matrix_core_dirichlet.tpl')
        with open(filename, 'r') as f:
            src = Template(f.read())
            code = src.substitute({
                'name': name.lower(),
                'init': 'nosh::matrix_core_dirichlet(%s)' % init,
                'eval_return_value': extract_c_expression(result),
                'eval_body': '' if is_x_used_eval else '(void) x;\n',
                })

        return {
            'code': code,
            }
