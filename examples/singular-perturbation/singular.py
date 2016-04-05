# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class Alpha(Expression):
    @classmethod
    def eval(cls, x): return sin(x[1])
    degree = 0


class Beta(Expression):
    pass


class Singular(FvmMatrix2):
    def eval(u):
        return \
              integrate(lambda x: n_dot(Beta(x)), dS())  # + \
              # integrate(lambda x: n_dot_grad(u(x)), dS())  # + \
              # integrate(lambda x: dot(n, grad(u(x))), dS()) #+ \
              # integrate(lambda x: u(x), dV())
              # integrate(lambda x: 2 * sin(x[0]) * u(x), dV())
    boundary_conditions = [Bc1()]


# class Core0(MatrixCore):
#     def edge_contrib(self, x0, x1, edge_length, edge_covolume):
#         eps = 2.0e-1
#         alpha = edge_covolume / edge_length
#         return [[eps * alpha, -eps * alpha],
#                 [-eps * alpha, eps * alpha]
#                 ]
#
#     def vertex_contrib(self, x, control_volume):
#         return control_volume
#
# class Singular(FvmMatrix):
#     matrix_cores = [Core0()]
#     boundary_conditions = [Bc1()]
