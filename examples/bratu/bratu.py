# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class NLaplace(FvmMatrix):
    def apply(u):
        return integrate(lambda x: -n_dot_grad(u, x), dS)

alpha = 0.75
# alpha = ScalarParameter()


class F(FvmOperator):
    def apply(u):
        return NLaplace(u) \
            - integrate(lambda x: alpha * exp(u(x)), dV)

    dirichlet = [(lambda x, u: u, Boundary)]


# class Jacobian(FvmOperator):
#     def apply(u, u0):
#         return NLaplace(u) \
#             - integrate(lambda x: alpha * exp(u0(x)) * u(x), dV)
#             # - integrate(lambda x: beta, dV)
#
#     dirichlet = [(lambda x, u: u, Boundary)]


# class Preconditioner(VectorOperator):
#     def eval(U, U0): return - Laplace(U)
#     def boundary_eval(U): return U


# class dFdp(FvmOperator):
#     def apply(u): return - integrate(lambda x: exp(u(x)), dV())
#     dirichlet = [(lambda x, u: 0.0, Boundary())]


# class Bratu(NonlinearFvmProblem):
#     F = F()
