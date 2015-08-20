# -*- coding: utf-8 -*-
# @HEADER
#
#    Operators for the Poisson problem.
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
import sympy

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


class A(EdgeMatrix):
    def edge_function(alpha, c0, c1):
        return [[alpha, -alpha], [-alpha, alpha]]
