# -*- coding: utf-8 -*-
from nfl import *
from sympy import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class F(Expression):
    def eval(x): return sin(x[1])
    degree = 0


class Laplace(FvmMatrix):
    def eval(u):
        return dot(n, grad(u)) * dS + dot(n, a) * u * dS

    boundary_conditions = [Bc1()]
