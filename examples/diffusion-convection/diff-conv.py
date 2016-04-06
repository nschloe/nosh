# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class F(Expression):
    def eval(x): return sin(x[1])
    degree = 0


class DC(FvmMatrix2):
    def eval(u):
        return integrate(
            lambda x: -n_dot_grad(u, x) + dot(n, Matrix([-1, -1, 0])) * u(x),
            dS()
        )
    boundary_conditions = [Bc1()]
