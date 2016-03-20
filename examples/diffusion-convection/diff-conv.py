# -*- coding: utf-8 -*-
from nfl import *
from sympy import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class F(Expression):
    def eval(x): return sin(x[1])
    degree = 0


class Laplace(FvmMatrix):
    def edge_contrib(x0, x1, edge_length, edge_covolume):
        a = [1, 0, 0]
        edge_midpoint = 0.5 * (x0 + x1)
        n = (x1 - x0) / edge_length
        beta = edge_covolume / 2 * (n.T * a)
        return [
            [
                alpha + beta,
                -alpha + beta
            ],
            [
                -alpha - beta,
                alpha - beta
            ]
            ]

    boundary_conditions = [Bc1()]
