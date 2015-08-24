# -*- coding: utf-8 -*-
from nfl import *
from sympy import *

bc1 = DirichletBC(
  lambda x: x[0] < 0,
  lambda x: 0.0
  )
bc2 = DirichletBC(
  lambda x: x[0] >= 0,
  lambda x: 1.0
  )

f = Expression(
    lambda x: sin(x[1]),
    degree=0
    )

# Laplace
class Laplace(FvmMatrix):
    def edge_contrib(alpha, edge_midpoint):
        return [[alpha, -alpha], [-alpha, alpha]]
