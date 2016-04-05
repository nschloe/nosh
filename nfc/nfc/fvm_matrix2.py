# -*- coding: utf-8 -*-
#
import nfl
import os
import sympy
from .helpers import extract_c_expression
from .discretize_edge_integral import DiscretizeEdgeIntegral


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
        print(res.vertex)
        v = self.get_code_vertex(u, res.vertex)
        print(v)
        e = self.get_code_edge(u, res.edge)
        print(e)
        return

    def get_code_vertex(self, u, function):
        # Numerically integrate function over a control volume.
        # TODO find out if the code can be evaluated at the volume "midpoint",
        # then evaluate it there, and multiply by the volume.
        x = sympy.MatrixSymbol('x', 1, 3)
        # Evaluate the function for u at x.
        print(function, u.is_Function)
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
        print('fu0', fu0)
        assert(is_linear(fu0, [u0]))
        coeff = fu0 / u0
        control_volume = sympy.Symbol('control_volume')
        coeff = control_volume * coeff
        print(coeff)
        return extract_c_expression(coeff)

    def get_code_edge(self, u, function):
        x = sympy.MatrixSymbol('x', 1, 3)

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
        coeff0 = expr.subs(generator.u1, 0) / generator.u0
        if generator.u0 in coeff0.free_symbols:
            coeff0 = sympy.simplify(coeff0)

        coeff1 = expr.subs(generator.u0, 0) / generator.u1
        if generator.u1 in coeff1.free_symbols:
            coeff1 = sympy.simplify(coeff1)

        exit()

        return ''
