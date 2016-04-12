# -*- coding: utf-8 -*-
from nfl import *
from sympy import sin


class Gamma0(Subdomain):
    def is_inside(self, x): return x[1] < 0
    is_boundary_only = True


class Gamma1(Subdomain):
    def is_inside(self, x): return x[1] >= 0
    is_boundary_only = True


class Poisson(LinearFvmProblem):
    def apply(u):
        return integrate(lambda x: -n_dot_grad(u, x), dS) \
                - integrate(lambda x: sin(x[1]), dV)

    dirichlet = [
            (lambda x: 0.0, Gamma0),
            (lambda x: 1.0, Gamma1)
            ]


# Alternative (raw) syntax:
# class Core0(MatrixCore):
#     def edge_contrib(self, x0, x1, edge_length, edge_covolume):
#         alpha = edge_covolume / edge_length
#         return [[alpha, -alpha], [-alpha, alpha]]
# class Laplace(FvmMatrix):
#     matrix_cores = [Core0()]
#     boundary_conditions = [Bc1(), Bc2()]
# class F(Expression):
#     def eval(x): return sin(x[1])
#     degree = 0
