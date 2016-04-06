# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class DC(LinearFvmProblem):
    def eval(u):
        a = Matrix([-1, -1, 0])
        return \
            integrate(
                lambda x: -n_dot_grad(u, x) + dot(n, a) * u(x),
                dS()
            ) \
            - integrate(lambda x: 1.0, dV())
    dirichlet_boundary_conditions = [Bc1()]
