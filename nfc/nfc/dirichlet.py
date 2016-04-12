# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy

from .helpers import \
        compare_variables, \
        extract_c_expression, \
        get_uuid, \
        sanitize_identifier, \
        templates_dir
from .subdomain import *


class Dirichlet(object):
    def __init__(self, function, subdomains, is_matrix):
        self.class_name = 'dirichlet_' + get_uuid()
        self.function = function
        self.dependencies = [SubdomainCode(sd) for sd in subdomains]
        if subdomains:
            self.subdomains = subdomains
        else:
            self.subdomains = [Boundary]
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_object(self, dep_class_objects):
        x = sympy.MatrixSymbol('x', 3, 1)
        result = self.function(x)

        unused_args, _ = compare_variables(set([x]), [result])

        # collect subdomain init code
        init = '{%s}' % ', '.join(
                '"%s"' % sanitize_identifier(sd.__name__)
                for sd in self.subdomains
                )

        # template substitution
        filename = os.path.join(templates_dir, 'matrix_core_dirichlet.tpl')
        with open(filename, 'r') as f:
            code = Template(f.read()).substitute({
                'name': self.class_name,
                'init': 'nosh::matrix_core_dirichlet(%s)' % init,
                'eval_return_value': extract_c_expression(result),
                'eval_body':
                    '\n'.join('(void) %s;' % arg for arg in unused_args)
                })

        return {
            'code': code,
            }
