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
import numpy

class Function(object):
    def __init__(self):
        return

class Coefficient(object):
    def __init__(self):
        return

class Expression(object):
    def __init__(self, expr, degree=numpy.infty):
        self.expr = expr
        return

class Vector(object):
    def __init__(self):
        return

class OuterNormal(Vector):
    def __init__(self):
        return

# Linear operator deifned via the operation along edges
# in a Delaunay mesh
class EdgeOperator(object):
    def __init__(self, edgeFunction):
        self.edgeFunction = edgeFunction
        #print(edgeFunction())
        return

class EdgeCoefficient(object):
    def __init__(self):
        return

class Integral(object):
    def __init__(self, integrand, measure):
        self.integrand = integrand
        self.measure = measure
        return

class DirichletBC(object):
    def __init__(self, insideCondition, evalReturn):
        self.isInside = insideCondition
        self.eval = evalReturn
        return

def inner(a, b):
    assert(isinstance(a, Vector))
    assert(isinstance(b, Vector))
    return

def dot(a, b):
    return inner(a, b)

def grad(a):
    return Vector()
