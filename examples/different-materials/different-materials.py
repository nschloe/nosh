# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class eps(Expression):
    pass


class Problem(LinearFvmProblem):
    def eval(u):
        return integrate(lambda x: -eps(x) * n_dot_grad(u, x), dS()) \
                + integrate(lambda x: -1.0, dV())

    dirichlet_boundary_conditions = [Bc1()]


# Alternative (raw) syntax:
# class Core(MatrixCore):
#     def edge_contrib(self, x0, x1, edge_length, edge_covolume):
#         alpha = edge_covolume / edge_length
#         edge_midpoint = 0.5 * (x0 + x1)
#         return [
#                 [
#                     eps(edge_midpoint) * alpha,
#                     -eps(edge_midpoint) * alpha
#                 ],
#                 [
#                     -eps(edge_midpoint) * alpha,
#                     eps(edge_midpoint) * alpha
#                 ]
#                 ]
#
#
# class Laplace(FvmMatrix):
#     matrix_cores = [Core()]
#     boundary_conditions = [Bc1()]
