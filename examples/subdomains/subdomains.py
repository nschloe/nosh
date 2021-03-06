# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class D1(Subdomain):
    def is_inside(self, x):
        return x[0]**2 + x[1]**2 > 25


class D2(Subdomain):
    def is_inside(self, x):
        return x[0]**2 + x[1]**2 <= 25


class P(LinearFvmProblem):
    def apply(u):
        return integrate(lambda x: -n_dot_grad(u(x)), dS, D1) \
            + integrate(lambda x: -0.2 * n_dot_grad(u(x)), dS, D2) \
            + integrate(lambda x: u(x), dV, D2) \
            - integrate(lambda x: 1.0, dV, D2)
    dirichlet = [(lambda x: 0.0, Boundary)]
