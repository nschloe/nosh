# -*- coding: utf-8 -*-
from nfl import *
from sympy import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


# class Singular(Operator):
#     def eval(u):
#         return - 2.0e-1 * dot(n, grad(u)) * domega + u * dV
#
#     boundary_conditions = [Bc1()]


class Singular(FvmMatrix):
    def edge_contrib(x0, x1, edge_length, edge_covolume):
        eps = 2.0e-1
        alpha = edge_covolume / edge_length
        return [[eps * alpha, -eps * alpha],
                [-eps * alpha, eps * alpha]
                ]

    def vertex_contrib(control_volume):
        return control_volume

    boundary_conditions = [Bc1()]
