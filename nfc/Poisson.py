# -*- coding: utf-8 -*-
# @HEADER
#
#    <description>
#    Copyright (C) 2015  Nico Schlömer
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# @HEADER
#
from nfl import *

#u = Function()
#f = Coefficient()
#f = Expression('x[0]*x[0] + 2*x[0] + 1', degree=2)

#n = OuterNormal()

#a = Integral(dot(n, grad(u)), 'ds')

bc1 = DirichletBC(
  lambda x: x[0] < 0,
  lambda x: 0.0
  )
bc2 = DirichletBC(
  lambda x: x[0] >= 0,
  lambda x: 1.0
  )
bc3 = DirichletBC(
  lambda x: x[0] > -10000,
  lambda x: 0.0
  )

f = Expression(
    lambda x: x[0]**2,
    degree=2
    )

def laplace(alpha, c0, c1):
  import sympy
  return [[alpha, -alpha],
          [-alpha, alpha]]
a = EdgeOperator(laplace)

#def glEdge(edge, alpha):
#  return [[alpha, -exp(1j * Integral(dot(edge, A), edge) * alpha],
#          [-exp(-1j * Integral(dot(edge, A), edge) * alpha, alpha]]
#a = EdgeOperator(glEdge)

#L = Integral(f, 'dx')