# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class D1(Subdomain):
    def is_inside(self, x):
        return x[0] < 0.0
    is_boundary_only = True


class Problem(LinearFvmProblem):
    def eval(u):
        return integrate(lambda x: -n_dot_grad(u, x), dS()) \
                + integrate(lambda x: 3.0, dGamma()) \
                - integrate(lambda x: 1.0, dV())

    dirichlet = [(lambda x: 0.0, D1())]
