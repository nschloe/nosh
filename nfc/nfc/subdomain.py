# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
from .helpers import extract_c_expression, get_uuid, templates_dir

import nfl


class SubdomainCode(object):
    def __init__(self, cls):
        self.cls = cls
        self.class_name = cls.__name__
        return

    def get_dependencies(self):
        return set()

    def get_class_obj(self):
        if self.cls == nfl.Boundary:
            # 'Boundary' is already defined
            return '', set()

        x = sympy.MatrixSymbol('x', 3, 1)

        result = obj.is_inside(x)

        expr_arguments = set([x])
        try:
            used_vars = result.free_variables
        except AttributeError:
            used_vars = set()
        unused_arguments = expr_arguments - used_vars

        # No undefined variables allowed
        assert(len(used_vars - expr_arguments) == 0)

        try:
            ibo = 'true' if subdomain.is_boundary_only else 'false'
        except AttributeError:
            # AttributeError: 'D2' object has no attribute 'is_boundary_only'
            ibo = 'false'

        # template substitution
        with open(os.path.join(templates_dir, 'subdomain.tpl'), 'r') as f:
            src = Template(f.read())
            code = src.substitute({
                'name': name.lower(),
                'id': '"%s"' % name.lower(),
                'is_inside_return': extract_c_expression(result),
                'is_boundary_only': ibo,
                'is_inside_body': '\n'.join(
                    ('(void) %s;' % name) for name in unused_arguments
                    ),
                })

        return {
            'code': code
            }
