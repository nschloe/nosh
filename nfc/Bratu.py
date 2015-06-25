# -*- coding: utf-8 -*-
# @HEADER
#
#    NFL file for the Bratu problem.
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

bc0 = DirichletBC(
    lambda x: x[0] > -10000,
    lambda x: 0.0
    )

laplace = EdgeMatrix(
    lambda alpha, c0, c1: [[alpha, -alpha],  [-alpha, alpha]]
    )

def f(u):
    return laplace(u) - lmbd * exp(u)
bratu = NonlinearOperator(
    evalu=f
    )
# Jac=(lambda u0: A - lmbd * exp(u0))
