# -*- coding: utf-8 -*-
#
import nfl
import os
import sympy
from .helpers import extract_c_expression
from .discretize_edge_integral import DiscretizeEdgeIntegral


class CodeFvmMatrix2(object):
    def __init__(self, obj, name):
        self.obj = obj
        self.name = name
        return

    def get_dependencies(self):
        return []

    def is_linear(self, expr, var):
        # Check if an expression is linear in a given var.
        if not sympy.Eq(sympy.diff(expr, var, var), 0):
            return False
        else:
            return True

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
        assert(self.is_linear(fu0, u0))
        coeff = fu0 / u0
        control_volume = sympy.Symbol('control_volume')
        coeff = control_volume * coeff
        print(coeff)
        return extract_c_expression(coeff)

    def get_code_edge(self, u, function):
        x = sympy.MatrixSymbol('x', 1, 3)

        generator = DiscretizeEdgeIntegral()
        expr = generator.generate(function(x))

        print(expr)

        exit()

        # # Evaluate the function for u at x.
        # fx = function(x)
        # print(fx)
        # exit()

        # # First, get rid of all differential operators; approximate them by
        # # finite differences.
        # u0 = sympy.Symbol('u0')
        # u1 = sympy.Symbol('u1')
        # # [...] TODO

        # if u in res.free_symbols:
        #     # Replace u by the value at the midpoint, 0.5 * (u0 + u1).
        #     pass

        # coedge_length = sympy.Symbol('coedge_length')
        # expr = coedge_length * 0.5 * (u0 + u1)

        # exit()
        return ''
