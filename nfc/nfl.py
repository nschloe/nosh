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
import sympy


class Function(sympy.Symbol):
    def __init__(self):
        return


class Coefficient(object):
    def __init__(self):
        return


class Expression(object):
    def __init__(self, eval, degree=sympy.oo):
        self.eval = eval
        self.degree = degree
        return


class Vector(object):
    def __init__(self):
        return


class ScalarParameter(object):
    def __init__(self, name):
        self.name = name
        return


class OuterNormal(Vector):
    def __init__(self):
        return


# Linear operator deifned via the operation along edges in a Delaunay mesh
class EdgeMatrix(sympy.Function):
    pass
    #def __call__(self, *args):
    #    print(super(EdgeMatrix, self))
    #    return super(EdgeMatrix, self).__call__(self, *args)
    #    print('Called !')

    ## Override __new__,
    ## cf. <https://groups.google.com/forum/#!topic/sympy/5mLEq4Gbyfk>.
    #def __new__(self, name, edge_function):
    #    obj = super(EdgeMatrix, self).__new__(self, name)
    #    return obj

    #def __init__(self, name, edge_function):
    #    super(EdgeMatrix, self).__init__(name)
    #    self.edge_function = edge_function
    #    return

    #def __new__(self, *args, **kwargs):
    #    #obj = super(EdgeMatrix, self).__new__(self, name)
    #    obj = super(EdgeMatrix, self).__new__(self, *args, **kwargs)
    #    return obj

    #def __init__(self, name):
    #    #obj = super(EdgeMatrix, self).__new__(self, name)
    #    #super(EdgeMatrix, self).__init__(self, name)
    #    return

    #def __init__(self, name, edge_function):
    #    super(EdgeMatrix, self).__init__(name)
    #    self.edge_function = edge_function
    #    return


class MatrixFactory(object):
    def __init__(self):
        return


class NonlinearOperator(object):
    def __init__(self, f=None, dfdp=None, jac=None):
        self.f = f
        self.dfdp = dfdp
        self.jac = jac
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
