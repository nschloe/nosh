# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .discretize_edge_integral import DiscretizeEdgeIntegral
from .code_generator_eigen import CodeGeneratorEigen
from .helpers import extract_c_expression, templates_dir


def _is_affine_linear(expr, vars):
    for var in vars:
        if not sympy.Eq(sympy.diff(expr, var, var), 0):
            return False
    return True


# We still need this for pure matrices
# def _is_linear(expr, vars):
#     if not _is_affine_linear(expr, vars):
#         return False
#     # Check that expr is not affine.
#     if isinstance(expr, int) or isinstance(expr, float):
#         return expr == 0
#     else:
#         return expr.subs([(var, 0) for var in vars]) == 0


def _get_expr_edge_contrib(self, method):
    x0 = sympy.MatrixSymbol('x0', 3, 1)
    x1 = sympy.MatrixSymbol('x1', 3, 1)
    edge_length = sympy.Symbol('edge_length')
    edge_covolume = sympy.Symbol('edge_covolume')
    all_symbols = set([x0, x1, edge_length, edge_covolume])

    specs = inspect.getargspec(method)
    assert(len(specs.args) == len(all_symbols) + 1)

    result = method(x0, x1, edge_length, edge_covolume)

    assert 'coeff' in result
    assert(len(result['coeff']) == 2)
    assert(len(result['coeff'][0]) == 2)
    assert(len(result['coeff'][1]) == 2)
    assert 'affine' in result
    assert(len(result['affine']) == 2)

    return result['coeff'], result['affine']


def get_core_code_generator(namespace, name, core):
    '''Get code generator from raw core object.
    '''
    # handle the edge contributions
    if callable(getattr(core, 'edge_contrib', None)):
        edge_coeff, edge_affine = \
                _get_expr_edge_contrib(core.edge_contrib)

    return CodeEdgeCore(
            namespace, name,
            edge_coeff, edge_affine
            )


def get_edge_core_code_from_expression(namespace, name, u, expr):
    '''Get code generator from expression.
    '''
    edge_coeff, edge_affine = _get_expr_edge(u, expr)

    return CodeEdgeCore(
            namespace, name,
            edge_coeff, edge_affine
            )


def _get_expr_edge(u, function):
    x = sympy.MatrixSymbol('x', 3, 1)

    generator = DiscretizeEdgeIntegral()
    expr, _ = generator.generate(function(x), u, x)

    # Make sure the expression is linear in u0, u1.
    if not _is_affine_linear(expr, [generator.u0, generator.u1]):
        raise RuntimeError((
            'The given function\n'
            '    f(x) = %s\n'
            'does not seem to be affine linear in u.')
            % function(x)
            )

    # Get the coefficients of u0, u1.
    coeff00 = sympy.diff(expr, generator.u0)
    coeff01 = sympy.diff(expr, generator.u1)

    # Now construct the coefficients for the other way around.
    coeff10 = coeff01.subs([
        (generator.u0, generator.u1),
        (generator.u1, generator.u0),
        (nfl.n, nfl.neg_n)
        ])
    coeff11 = coeff00.subs([
        (generator.u0, generator.u1),
        (generator.u1, generator.u0),
        (nfl.n, nfl.neg_n)
        ])

    affine = expr.subs([(generator.u0, 0), (generator.u1, 0)])

    return (
        [[coeff00, coeff01], [coeff10, coeff11]],
        [affine, affine]
        )


class CodeEdgeCore(object):
    def __init__(
            self, namespace, name,
            edge_coeff, edge_affine
            ):
        # edge
        edge_body, used_expressions = \
            self.get_code_elements_edge(edge_coeff, edge_affine)

        # now take care of the template substitution
        members_init = []
        members_declare = []

        self.dependencies = set()
        self.dependencies.update(used_expressions)

        # init and declare expressions in the C++ code
        for expr in used_expressions:
            members_init.append('%s(%s::%s())' % (expr, self.namespace, expr))
            members_declare.append(
                    'const %s::%s %s;' % (self.namespace, expr, expr)
                    )

        if members_init:
            members_init_code = ':\n' + ',\n'.join(members_init)
        else:
            members_init_code = ''

        self.class_name = name.lower() + '_edge_core'
        # template substitution
        with open(os.path.join(templates_dir, 'edge_core.tpl'), 'r') as f:
            src = Template(f.read())
            self.code = src.substitute({
                'name': self.class_name,
                'edge00': self.expr_to_code(edge_coeff[0][0]),
                'edge01': self.expr_to_code(edge_coeff[0][1]),
                'edge10': self.expr_to_code(edge_coeff[1][0]),
                'edge11': self.expr_to_code(edge_coeff[1][1]),
                'edge_affine0': self.expr_to_code(-edge_affine[0]),
                'edge_affine1': self.expr_to_code(-edge_affine[1]),
                'edge_body': '\n'.join(edge_body),
                'members_init': members_init_code,
                'members_declare': '\n'.join(members_declare)
                })
        return

    def get_dependencies(self):
        return self.dependencies

    def get_code(self):
        return self.code

    def expr_to_code(self, expr):
        gen = CodeGeneratorEigen()
        code, _ = gen.generate(expr)
        return code

    def get_code_elements_edge(self, edge_coeff, edge_affine):
        # Check out used an unused symbols for the edge code.
        # Arguments from the template.
        arguments = set([
            sympy.Symbol('x0'),
            sympy.Symbol('x1'),
            sympy.Symbol('edge_length'),
            sympy.Symbol('edge_covolume')
            ])
        expressions = [
            edge_coeff[0][0], edge_coeff[0][1],
            edge_coeff[1][0], edge_coeff[1][1],
            edge_affine[0], edge_affine[1]
            ]

        used_symbols = set([])
        used_expressions = set([])
        body = []
        for expr in expressions:
            try:
                used_symbols.update(expr.free_symbols)
            except AttributeError:
                pass

            used_expressions.update(set([
                    type(atom) for atom in expr.atoms(nfl.Expression)
                    ]))

        unused_arguments = arguments - used_symbols
        undefined_symbols = used_symbols - arguments

        # special treatment for n
        if nfl.n in undefined_symbols or \
           nfl.neg_n in undefined_symbols:
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

        return body, used_expressions

    def scan_code(self, arguments, expressions):
        used_symbols = set([])
        used_expressions = set([])

        for expr in expressions:
            try:
                used_symbols.update(expr.free_symbols)
            except AttributeError:
                pass

            used_expressions.update(set([
                    type(atom) for atom in expr.atoms(nfl.Expression)
                    ]))

        unused_arguments = arguments - used_symbols
        undefined_symbols = used_symbols - arguments
        assert(len(undefined_symbols) == 0)

        return unused_arguments, used_expressions
