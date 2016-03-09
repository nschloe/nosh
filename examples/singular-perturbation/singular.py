# -*- coding: utf-8 -*-
from nfl import *
from sympy import *


class Bc1(DirichletBC):
    def eval(self, x): return 0.0


class Singular(FvmMatrix):
    def edge_contrib(alpha, edge_midpoint):
        eps = 2.0e-1
        return [[eps * alpha, -eps * alpha],
                [-eps * alpha, eps * alpha]
                ]

    def vertex_contrib(control_volume):
        return control_volume

    boundary_conditions = [Bc1()]
