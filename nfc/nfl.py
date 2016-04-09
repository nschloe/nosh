# -*- coding: utf-8 -*-
#
import sympy


class Coefficient(object):
    pass


class Operator(sympy.Function):
    pass


class Expression(sympy.Function):
    pass
    # degree = sympy.oo


class Vector(object):
    pass


class ScalarParameter(object):
    pass


class OuterNormal(Vector):
    pass


class Subdomain(object):
    pass


class FvmMatrix(sympy.Function):
    pass


class LinearFvmProblem(object):
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


class dGamma(Measure):
    pass


def integrate(integrand, measure, subdomains=None):
    assert(isinstance(measure, Measure))

    if subdomains is None:
        subdomains = set()
    elif not isinstance(subdomains, set):
        try:
            subdomains = set(subdomains)
        except TypeError:  # TypeError: 'D1' object is not iterable
            subdomains = set([subdomains])

    assert(
        isinstance(measure, dS) or
        isinstance(measure, dV) or
        isinstance(measure, dGamma)
        )

    return Core([(integrand, measure, subdomains)])


class Core(object):
    def __init__(self, set_c):
        self.cores = set_c

    def __add__(self, other):
        self.cores.extend(other.cores)
        return Core(self.cores)

    # def __sub__(self, other):
    #     return Core(
    #             lambda x: self.vertex(x) - other.vertex(x),
    #             lambda x: self.edge(x) - other.edge(x),
    #             lambda x: self.domain_boundary(x) - other.domain_boundary(x)
    #             )


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
    def eval(self, x): return 0.0
    subdomains = set()


class NeumannBC(object):
    def is_inside(self, x): return x[0] < 1e10

    def eval(self, x): return 0.0


class EdgeCore(object):
    pass


def inner(a, b):
    assert(isinstance(a, Vector))
    assert(isinstance(b, Vector))
    return


class dot(sympy.Function):
    pass


class n_dot_grad(sympy.Function):
    pass

n = sympy.MatrixSymbol('n', 3, 1)
neg_n = sympy.MatrixSymbol('neg_n', 3, 1)


def grad(a):
    return Vector()
