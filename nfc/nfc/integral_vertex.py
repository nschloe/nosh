# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .code_generator_eigen import CodeGeneratorEigen
from .expression import *
from .subdomain import *
from .helpers import \
        extract_c_expression, \
        extract_linear_components, \
        get_uuid, \
        is_affine_linear, \
        members_init_declare, \
        templates_dir


class IntegralVertex(object):
    def __init__(self, namespace, u, integrand, subdomains, is_matrix):
        self.namespace = namespace
        self.class_name = 'vertex_core_' + get_uuid()

        self.is_matrix = is_matrix

        self.expr, self.u0 = \
            _discretize_integral(u, integrand)

        self.dependencies = set().union(
            [ExpressionCode(type(atom))
                for atom in self.expr.atoms(nfl.Expression)],
            [SubdomainCode(sd) for sd in subdomains]
            )
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_object(self, dependency_class_objects):
        arguments = set([
            sympy.MatrixSymbol('x', 3, 1),
            sympy.Symbol('control_volume')
            ])
        used_vars = self.expr.free_symbols
        if self.u0 in used_vars:
            used_vars.remove(self.u0)
        unused_args = arguments - used_vars

        if self.is_matrix:
            parent_class = 'matrix_core_vertex'
        else:
            parent_class = 'operator_core_vertex'

        # now take care of the template substitution
        members_init, members_declare = \
            members_init_declare(
                    self.namespace,
                    parent_class,
                    dependency_class_objects
                    )

        if members_init:
            members_init_code = ':\n' + ',\n'.join(members_init)
        else:
            members_init_code = ''

        if self.is_matrix:
            coeff, affine = extract_linear_components(self.expr, self.u0)
            type = 'matrix_core_vertex'
            filename = os.path.join(templates_dir, 'matrix_core_vertex.tpl')
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'vertex_contrib': extract_c_expression(coeff),
                    'vertex_affine': extract_c_expression(-affine),
                    'vertex_body': '\n'.join(
                        ('(void) %s;' % name) for name in unused_args
                        ),
                    'members_init': members_init_code,
                    'members_declare': '\n'.join(members_declare)
                    })
        else:
            type = 'operator_core_vertex'
            filename = os.path.join(templates_dir, 'operator_core_vertex.tpl')
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'return_value': extract_c_expression(self.expr),
                    'eval_body': '\n'.join(
                        ('(void) %s;' % name) for name in unused_args
                        ),
                    'members_init': members_init_code,
                    'members_declare': '\n'.join(members_declare)
                    })

        return {
            'type': type,
            'code': code,
            'class_name': self.class_name,
            'constructor_args': []
            }


# def get_matrix_core_vertex_code(namespace, class_name, core):
#     '''Get code generator from raw core object.
#     '''
#     # handle the vertex contributions
#     x = sympy.MatrixSymbol('x')
#     vol = sympy.Symbol('control_volume')
#     all_symbols = set([x, vol])
#
#     specs = inspect.getargspec(method)
#     assert(len(specs.args) == len(all_symbols) + 1)
#
#     vertex_coeff, vertex_affine = method(x, vol)
#
#     return _get_code_matrix_core_vertex(
#             namespace, class_name,
#             vertex_coeff, vertex_affine
#             )


def _discretize_integral(u, function):
    # Numerically integrate function over a control volume.
    x = sympy.MatrixSymbol('x', 3, 1)
    # Evaluate the function for u at x.
    fx = function(x)
    # Replace all occurences of u(x) by u0 (the value at the control volume
    # center) and multiply by the control volume)
    u0 = sympy.Symbol('u0')
    try:
        fu0 = fx.subs(u(x), u0)
    except AttributeError:
        # 'int' object has no attribute 'subs'
        fu0 = fx
    control_volume = sympy.Symbol('control_volume')
    return control_volume * fu0, u0
