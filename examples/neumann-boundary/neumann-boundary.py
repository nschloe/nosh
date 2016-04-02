# -*- coding: utf-8 -*-
from nfl import *
from sympy import *


class Bc1(DirichletBC):
    def is_inside(self, x): return x[0] < 0 and x[1] < 0

    def eval(self, x): return 0.0


class Bc2(DirichletBC):
    def is_inside(self, x): return x[0] < 0 and x[1] >= 0

    def eval(self, x): return 1.0


class Bc3(NeumannBC):
    def is_inside(self, x): return x[0] >= 0

    # Careful! This must be seen in context with edge_contrib of whatever
    # operator it's used with.
    # Here: Neumann condition with n.grad(u) = 1.
    def eval(self, x, surface): return 1.0 * surface


class F(Expression):
    def eval(x): return 1.0
    degree = 0


class Laplace(FvmMatrix):
    def edge_contrib(x0, x1, edge_length, edge_covolume):
        alpha = edge_covolume / edge_length
        return [[alpha, -alpha], [-alpha, alpha]]

    boundary_conditions = [Bc1(), Bc2(), Bc3()]

