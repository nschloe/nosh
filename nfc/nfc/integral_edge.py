# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .discretize_edge_integral import discretize_edge_integral
from .code_generator_eigen import CodeGeneratorEigen
from .helpers import \
        extract_c_expression, \
        is_affine_linear, \
        get_uuid, \
        members_init_declare, \
        templates_dir


class IntegralEdge(object):
    def __init__(self, u, integrand, subdomains, is_matrix):
        self.class_name = 'matrix_edge_core_' + get_uuid()

        self.expr, self.u0, self.u1 = \
            discretize_edge_integral(integrand, u)

        self.is_matrix = is_matrix

        # gather expressions and subdomains as dependencies
        self.dependencies = set().union(
            [type(atom) for atom in self.expr.atoms(nfl.Expression)],
            subdomains
            )
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_object(self, dependency_class_objects):
        # Arguments from the template.
        arguments = set([
            sympy.Symbol('x0'),
            sympy.Symbol('x1'),
            sympy.Symbol('edge_length'),
            sympy.Symbol('edge_covolume')
            ])

        used_vars = self.expr.free_symbols
        used_vars.remove(self.u0)
        used_vars.remove(self.u1)
        eval_body = _get_code_body(arguments, used_vars)

        members_init, members_declare = \
            members_init_declare('matrix_core_edge', dependency_class_objects)

        if members_init:
            members_init_code = ':\n' + ',\n'.join(members_init)
        else:
            members_init_code = ''

        if self.is_matrix:
            type = 'matrix_core_edge'
            edge_coeff, edge_affine = \
                _extract_linear_components(self.expr, self.u0, self.u1)
            # template substitution
            filename = os.path.join(templates_dir, 'matrix_core_edge.tpl')
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'edge00': _expr_to_code(edge_coeff[0][0]),
                    'edge01': _expr_to_code(edge_coeff[0][1]),
                    'edge10': _expr_to_code(edge_coeff[1][0]),
                    'edge11': _expr_to_code(edge_coeff[1][1]),
                    'edge_affine0': _expr_to_code(-edge_affine[0]),
                    'edge_affine1': _expr_to_code(-edge_affine[1]),
                    'eval_body': '\n'.join(eval_body),
                    'members_init': members_init_code,
                    'members_declare': '\n'.join(members_declare)
                    })
        else:
            # TODO
            type = 'operator_core_edge'
            # template substitution
            filename = os.path.join(templates_dir, 'operator_core_edge.tpl')
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'edge00': _expr_to_code(edge_coeff[0][0]),
                    'edge01': _expr_to_code(edge_coeff[0][1]),
                    'edge10': _expr_to_code(edge_coeff[1][0]),
                    'edge11': _expr_to_code(edge_coeff[1][1]),
                    'edge_affine0': _expr_to_code(-edge_affine[0]),
                    'edge_affine1': _expr_to_code(-edge_affine[1]),
                    'eval_body': '\n'.join(eval_body),
                    'members_init': members_init_code,
                    'members_declare': '\n'.join(members_declare)
                    })

        return {
            'type': type,
            'code': code,
            'class_name': self.class_name,
            'constructor_args': []
            }


def _extract_linear_components(expr, u0, u1):
    if not is_affine_linear(expr, [u0, u1]):
        raise RuntimeError((
            'The given function\n'
            '    f(x) = %s\n'
            'does not seem to be affine linear in u.')
            % function(x)
            )

    # Get the coefficients of u0, u1.
    coeff00 = sympy.diff(expr, u0)
    coeff01 = sympy.diff(expr, u1)

    # Now construct the coefficients for the other way around.
    coeff10 = coeff01.subs([(u0, u1), (u1, u0), (nfl.n, nfl.neg_n)])
    coeff11 = coeff00.subs([(u0, u1), (u1, u0), (nfl.n, nfl.neg_n)])

    affine = expr.subs([(u0, 0), (u1, 0)])

    return (
        [[coeff00, coeff01], [coeff10, coeff11]],
        [affine, affine]
        )


def _expr_to_code(expr):
    gen = CodeGeneratorEigen()
    code, _ = gen.generate(expr)
    return code


def _get_code_body(arguments, used_variables):
    unused_arguments = arguments - used_variables
    undefined_symbols = used_variables - arguments

    body = []
    # special treatment for n
    if nfl.n in undefined_symbols or nfl.neg_n in undefined_symbols:
        body.append('const auto n = (x1 - x0) / edge_length;')
        unused_arguments -= set([
            sympy.Symbol('x0'),
            sympy.Symbol('x1')
            ])
        undefined_symbols.remove(nfl.n)
        undefined_symbols.remove(nfl.neg_n)

    assert(len(undefined_symbols) == 0)

    for name in unused_arguments:
        body.append('(void) %s;' % name)

    return body
