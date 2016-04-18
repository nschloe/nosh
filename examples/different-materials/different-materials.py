# -*- coding: utf-8 -*-
from nfl import *


class eps(Expression):
    eval_body = '''
    if (x[0] > 0.0) {
      return 3.0;
    } else {
      return 1.0;
    }
    '''


class Problem(LinearFvmProblem):
    def apply(u):
        return integrate(lambda x: -eps(x) * n_dot_grad(u(x)), dS) \
                - integrate(lambda x: 1.0, dV)

    dirichlet = [(lambda x: 0.0, Boundary)]


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
