# -*- coding: utf-8 -*-
#
import sympy


class Function(sympy.Symbol):
    def __init__(self):
        return


class Coefficient(object):
    def __init__(self):
        return


class Operator(sympy.Function):
    pass


class Expression(sympy.Function):
    pass
    # degree = sympy.oo


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


# Linear operator defined via the operation along edges in a Delaunay mesh
class FvmMatrix(sympy.Function):
    pass

    #def edge_contrib(self, alpha, edge_midpoint):
    #    return []

    #def vertex_contrib(self, control_volume):
    #    return None

    #def __call__(self, *args):
    #    print(super(FvmMatrix, self))
    #    return super(FvmMatrix, self).__call__(self, *args)
    #    print('Called !')

    ## Override __new__,
    ## cf. <https://groups.google.com/forum/#!topic/sympy/5mLEq4Gbyfk>.
    #def __new__(self, name, edge_function):
    #    obj = super(FvmMatrix, self).__new__(self, name)
    #    return obj

    #def __init__(self, name, edge_function):
    #    super(FvmMatrix, self).__init__(name)
    #    self.edge_function = edge_function
    #    return

    #def __new__(self, *args, **kwargs):
    #    #obj = super(FvmMatrix, self).__new__(self, name)
    #    obj = super(FvmMatrix, self).__new__(self, *args, **kwargs)
    #    return obj

    #def __init__(self, name):
    #    #obj = super(FvmMatrix, self).__new__(self, name)
    #    #super(FvmMatrix, self).__init__(self, name)
    #    return

    #def __init__(self, name, edge_function):
    #    super(FvmMatrix, self).__init__(name)
    #    self.edge_function = edge_function
    #    return

class FvmMatrix2(sympy.Function):
    pass


class MatrixFactory(object):
    def __init__(self):
        return


class Measure(object):
    pass


class dV(Measure):
    pass


class dS(Measure):
    pass


def integrate(integrand, measure):
    assert(isinstance(measure, Measure))

    if isinstance(measure, dV):
        return Contributions(integrand, 0)
    elif isinstance(measure, dS):
        return Contributions(0, integrand)
    else:
        raise RuntimeError('Illegal measure')


class Contributions(object):
    def __init__(self, vertex, edge):
        self.vertex = vertex
        self.edge = edge

    def __add__(self, other):
        return Contributions(
                lambda x: self.vertex(x) + other.vertex(x),
                lambda x: self.edge(x) + other.edge(x)
                )

    def __sub__(self, other):
        return Contributions(
                lambda x: self.vertex(x) - other.vertex(x),
                lambda x: self.edge(x) - other.edge(x)
                )


class NonlinearProblem(object):
    def __init__(self, f=None, dfdp=None, jac=None, prec=None):
        self.f = f
        self.dfdp = dfdp
        self.jac = jac
        self.prec = prec
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
    def is_inside(self, x): return x[0] < 1e10

    def eval(self, x): return 0.0


class NeumannBC(object):
    def is_inside(self, x): return x[0] < 1e10

    def eval(self, x): return 0.0


class MatrixCore(object):
    def is_inside(self, x): return x[0] < 1e10

    def edge_contrib(self, x0, x1, edge_length, edge_covolume):
        return [[0.0, 0.0], [0.0, 0.0]]

    def vertex_contrib(self, x, control_volume):
        return 0.0


def inner(a, b):
    assert(isinstance(a, Vector))
    assert(isinstance(b, Vector))
    return


def dot(a, b):
    return inner(a, b)


def grad(a):
    return Vector()
