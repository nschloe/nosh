# -*- coding: utf-8 -*-
from nfl import *
from sympy import *


class Bc1(DirichletBC):
    def is_inside(self, x): return x[1] < 0

    def eval(self, x): return 0.0


class Bc2(DirichletBC):
    def is_inside(self, x): return x[1] >= 0

    def eval(self, x): return 1.0


class F(Expression):
    def eval(x): return sin(x[1])
    degree = 0


class Laplace(FvmMatrix):
    def edge_contrib(x0, x1, edge_length, edge_covolume):
        alpha = edge_covolume / edge_length
        return [[alpha, -alpha], [-alpha, alpha]]

    boundary_conditions = [Bc1(), Bc2()]

# dS n grad(u)

# class Laplace(Operator):
#     def eval(u):
#         return integral(- dot(n, grad(u)), dS)
#
#     boundary_conditions = [Bc1(), Bc2()]
