# -*- coding: utf-8 -*-
from nfl import *
from sympy import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class F(Expression):
    def eval(self, x): return 1.0
    degree = 0


class Laplace(FvmMatrix):
    def edge_contrib(alpha, edge_midpoint):
        return [[alpha, -alpha], [-alpha, alpha]]
    boundary_conditions = [Bc1()]
