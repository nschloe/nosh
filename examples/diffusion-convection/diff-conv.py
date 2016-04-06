# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class DC(LinearFvmProblem):
    def eval(u):
        return \
            integrate(
                lambda x: -n_dot_grad(u, x) + dot(n, Matrix([-1, -1, 0])) * u(x),
                dS()
            ) \
            - integrate(lambda x: 1.0, dV())
    boundary_conditions = [Bc1()]
