# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class Bc1(DirichletBC):
    def is_inside(self, x): return x[0] < 0

    def eval(self, x): return 0.0


class Problem(LinearFvmProblem):
    def eval(u):
        return integrate(lambda x: -n_dot_grad(u, x), dS()) \
                + integrate(lambda x: 1.0, dGamma()) \
                - integrate(lambda x: 1.0, dV())
    dirichlet_boundary_conditions = [Bc1()]
