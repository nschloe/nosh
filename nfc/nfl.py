# -*- coding: utf-8 -*-
#
import sympy


class FvmOperator(sympy.Function):
    pass


class Expression(sympy.Function):
    pass
    # degree = sympy.oo


class Vector(object):
    pass


class ScalarParameter(object):
    pass


class Subdomain(object):
    pass


class Boundary(Subdomain):
    pass


class FvmMatrix(sympy.Function):
    # By default: No Dirichlet conditions.
    dirichlet = []


class LinearFvmProblem(object):
    pass


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

    return CoreList([Core(integrand, measure, subdomains)])


class Core(object):
    def __init__(self, integrand, measure, subdomains):
        self.integrand = integrand
        self.measure = measure
        self.subdomains = subdomains
        return


class CoreList(object):
    def __init__(self, list_cores):
        self.cores = list_cores

    def __add__(self, other):
        self.cores.extend(other.cores)
        return self

    def __sub__(self, other):
        # flip the sign on the integrand of all 'other' cores
        new_cores = []
        for core in other.cores:
            new_cores.append(Core(
                lambda x: -core.integrand(x),
                core.measure,
                core.subdomains
                ))
        self.cores.extend(new_cores)
        return self

    def __pos__(self):
        return self

    def __neg__(self):
        # flip the sign on the integrand of all 'self' cores
        new_cores = []
        for core in self.cores:
            new_cores.append(Core(
                lambda x: -core.integrand(x),
                core.measure,
                core.subdomains
                ))
        self.cores = new_cores
        return self

    def __mul__(self, other):
        assert(isinstance(other, float) or isinstance(other, int))
        # flip the sign on the integrand of all 'self' cores
        new_cores = []
        for core in self.cores:
            new_cores.append(Core(
                lambda x: other * core.integrand(x),
                core.measure,
                core.subdomains
                ))
        self.cores = new_cores
        return self

    __rmul__ = __mul__


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
