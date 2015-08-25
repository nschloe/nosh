# -*- coding: utf-8 -*-
from nfl import *
from sympy import *

bc1 = DirichletBC(
  lambda x: x[0] < 1e6,  # True
  lambda x: 0.0
  )


class Singular(FvmMatrix):
    def edge_contrib(alpha, edge_midpoint):
        eps = 2.0e-1
        return [[eps * alpha, -eps * alpha],
                [-eps * alpha, eps * alpha]
                ]

    def vertex_contrib(control_volume):
        return control_volume
