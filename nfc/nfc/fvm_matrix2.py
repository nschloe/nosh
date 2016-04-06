# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
import sympy
from .helpers import extract_c_expression, templates_dir
from .discretize_edge_integral import DiscretizeEdgeIntegral
from .code_generator_eigen import CodeGeneratorEigen


def is_affine_linear(expr, vars):
    for var in vars:
        if not sympy.Eq(sympy.diff(expr, var, var), 0):
            return False
    return True


def is_linear(expr, vars):
    if not is_affine_linear(expr, vars):
        return False
    # Check that expr is not affine.
    if isinstance(expr, int):
        return expr == 0
    else:
        return expr.subs([(var, 0) for var in vars]) == 0


class CodeFvmMatrix2(object):
    def __init__(self, obj, name):
        self.obj = obj
        self.name = name
        return

    def get_dependencies(self):
        return []

    def get_code(self):
        u = sympy.Function('u')
        res = self.obj.eval(u)
        assert(isinstance(res, nfl.Core))
        v = self.get_expr_vertex(u, res.vertex)
        e = self.get_expr_edge(u, res.edge)

        name = 'dummy'
        edge_unused_symbols = set([])
        vertex_unused_symbols = set([])
        members_init_code = ''
        members_declare = []

        # template substitution
        with open(os.path.join(templates_dir, 'matrix_core.tpl'), 'r') as f:
            src = Template(f.read())
            matrix_core_code = src.substitute({
                'name': name.lower(),  # class names are lowercase
                'edge00': self.expr_to_code(e[0][0]),
                'edge01': self.expr_to_code(e[0][1]),
                'edge10': self.expr_to_code(e[1][0]),
                'edge11': self.expr_to_code(e[1][1]),
                'edge_unused_args': '\n'.join(
                    ('(void) %s;' % name) for name in edge_unused_symbols
                    ),
                'vertex_contrib': extract_c_expression(v),
                'vertex_unused_args': '\n'.join(
                    ('(void) %s;' % name) for name in vertex_unused_symbols
                    ),
                'members_init': members_init_code,
                'members_declare': '\n'.join(members_declare)
                })
        print(matrix_core_code)
        exit()
        return

    def expr_to_code(self, expr):
        gen = CodeGeneratorEigen()
        code, _ = gen.generate(expr)
        return code

    def get_expr_vertex(self, u, function):
        # Numerically integrate function over a control volume.
        # TODO find out if the code can be evaluated at the volume "midpoint",
        # then evaluate it there, and multiply by the volume.
        x = sympy.MatrixSymbol('x', 1, 3)
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
        # Make sure that f is linear in u0; we're building a matrix here.
        assert(is_linear(fu0, [u0]))
        coeff = sympy.diff(fu0, u0)
        control_volume = sympy.Symbol('control_volume')
        coeff = control_volume * coeff
        return coeff

    def get_expr_edge(self, u, function):
        x = sympy.MatrixSymbol('x', 3, 1)

        generator = DiscretizeEdgeIntegral()
        expr, _ = generator.generate(function(x), u, x)

        # Make sure the expression is linear in u0, u1.
        if not is_linear(expr, [generator.u0, generator.u1]):
            raise RuntimeError((
                'The given function\n'
                '    f(x) = %s\n'
                'does not seem to be linear in u.')
                % function(x)
                )

        # Get the coefficients of u0, u1.
        coeff00 = sympy.diff(expr, generator.u0)
        coeff01 = sympy.diff(expr, generator.u1)

        # Now construct the coefficients for the other way around.
        coeff10 = coeff01.subs([
            (generator.u0, generator.u1),
            (generator.u1, generator.u0),
            (nfl.n, -nfl.n)
            ])
        coeff11 = coeff00.subs([
            (generator.u0, generator.u1),
            (generator.u1, generator.u0),
            (nfl.n, -nfl.n)
            ])

        return [[coeff00, coeff01], [coeff10, coeff11]]
