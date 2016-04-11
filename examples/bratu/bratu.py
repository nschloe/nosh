# -*- coding: utf-8 -*-
from sympy import *
from nfl import *


class NLaplace(FvmMatrix):
    def apply(u):
        return integrate(lambda x: -n_dot_grad(u, x), dS())

alpha = 0.75
# alpha = Parameter()


class F(FvmOperator):
    def apply(u):
        return NLaplace(u) - integrate(lambda x: alpha * exp(u(x)), dV())

    dirichlet = [(lambda x, u: u, Boundary())]


# class dFdp(Operator):
#     def eval(U): return -exp(U)
#     def boundary_eval(U): return 0.0
#
#
# class Jacobian(Operator):
#     U0 = Vector()
#     def eval(U): return - Laplace(U) - alpha * exp(U0) * U
#     def boundary_eval(U): return U


# class Preconditioner(VectorOperator):
#     def eval(U, U0): return - Laplace(U)
#     def boundary_eval(U): return U


# class Bratu(NonlinearFvmProblem):
#     F = F()



# class Bu(NonlinearDirichletBC):
#     def eval(u, x): return u(x)
#
# class B0(NonlinearDirichletBC):
#     def eval(u, x): return 0.0
#
# class F(VectorOperator):
#     def eval(U): return - Laplace(U) - alpha * exp(U)
#     boundary_evals = [Bu()]
#
# class dFdp(VectorOperator):
#     def eval(u): return - exp(u)
#     boundary_evals = [B0()]
#
# class Jacobian(VectorOperator):
#     def eval(u, u0): return - Laplace(u) - alpha * exp(u0) * u
#     boundary_evals = [Bu()]
#
# class Prec(VectorOperator):
#     def eval(u, u0): return - Laplace(u)
#     boundary_evals = [Bu()]
