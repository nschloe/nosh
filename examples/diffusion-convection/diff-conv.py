# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class DC(LinearFvmProblem):
    def eval(u):
        a = Matrix([-1, -1, 0])
        return \
            integrate(lambda x: -n_dot_grad(u, x) + dot(n, a) * u(x), dS()) \
            + integrate(lambda x: -1.0, dV())
    dirichlet = [(lambda x: 0.0, Boundary())]
