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
import sympy

bc0 = DirichletBC(
    lambda x: x[0] > -10000,
    lambda x: 0.0
    )


class laplace(EdgeMatrix):
    def edge_function(alpha, c0, c1):
        return [[alpha, -alpha], [-alpha, alpha]]

# lmbd = ScalarParameter('lambda')
lmbd = sympy.Symbol('lambda')
# alpha = 2.0
# beta = - 4.0
# f =  u + alpha + 1 + beta
# f =  u + 1 + 2 + sympy.exp(u / 6)
# f =  laplace(u) + 2
# f =  sympy.exp(u)
# f =  2 * u + laplace(u)
# f =  daplace(laplace(u))
# f = laplace(u) - lmbd * sympy.exp(u)

bratu = NonlinearOperator(
    f=lambda u: laplace(u) - lmbd * sympy.exp(u),
    dfdp=lambda u: -sympy.exp(u),
    jac=lambda u, u0: laplace(u) - lmbd * sympy.exp(u0) * u
    )
