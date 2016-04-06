# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class Bc1(DirichletBC):
    def is_inside(self, x): return x[1] < 0

    def eval(self, x): return 0.0


class Bc2(DirichletBC):
    def is_inside(self, x): return x[1] >= 0

    def eval(self, x): return 1.0


class F(Expression):
    def eval(x): return sin(x[1])
    degree = 0


class Laplace(FvmMatrix2):
    def eval(u):
        return integrate(
            lambda x: -n_dot_grad(u, x),
            dS()
        )
    boundary_conditions = [Bc1(), Bc2()]

# Alternative (raw) syntax:
# class Core0(MatrixCore):
#     def edge_contrib(self, x0, x1, edge_length, edge_covolume):
#         alpha = edge_covolume / edge_length
#         return [[alpha, -alpha], [-alpha, alpha]]
# class Laplace(FvmMatrix):
#     matrix_cores = [Core0()]
#     boundary_conditions = [Bc1(), Bc2()]
