# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class Singular(LinearFvmProblem):
    def eval(u):
        return integrate(lambda x: - 0.2 * n_dot_grad(u, x), dS()) \
               + integrate(lambda x: u(x), dV()) \
               + integrate(lambda x: -1.0, dV())
    dirichlet = [(lambda x: 0.0, Boundary())]


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
