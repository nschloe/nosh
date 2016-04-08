# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class D1(Subdomain):
    def is_inside(self, x):
        return x[1] < 0
    is_boundary_only = True


class D2(Subdomain):
    def is_inside(self, x):
        return x[1] >= 0
    is_boundary_only = True


class Bc1(DirichletBC):
    def eval(self, x): return 0.0
    subdomains = [D1()]


class Bc2(DirichletBC):
    def eval(self, x): return 1.0
    subdomains = [D2()]


class Poisson(LinearFvmProblem):
    def eval(u):
        return integrate(lambda x: -n_dot_grad(u, x), dS()) \
                + integrate(lambda x: - sin(x[1]), dV())
    dirichlet_boundary_conditions = [Bc1(), Bc2()]


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
