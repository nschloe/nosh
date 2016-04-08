# -*- coding: utf-8 -*-
#
import sympy


class Function(sympy.Symbol):
    def __init__(self):
        return


class Coefficient(object):
    pass


class VectorOperator(sympy.Function):
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


class LinearFvmProblem(object):
    pass


class MatrixFactory(object):
    def __init__(self):
        return


class Measure(object):
    def __init__(self, subdomains=None):
        if subdomains is None:
            self.subdomains = set()
        else:
            try:
                self.subdomains = set(subdomains)
            except TypeError:
                # TypeError: 'B' object is not iterable
                self.subdomains = set([subdomains])
        return


class dV(Measure):
    def __init__(self, subdomains=None):
        super().__init__(subdomains)


class dS(Measure):
    def __init__(self, subdomains=None):
        super().__init__(subdomains)


class dGamma(Measure):
    def __init__(self, subdomains=None):
        super().__init__(subdomains)


def integrate(integrand, measure):
    assert(isinstance(measure, Measure))

    if isinstance(measure, dS):
        return Core(set([(integrand, measure)]))
    elif isinstance(measure, dV):
        return Core(set([(integrand, measure)]))
    elif isinstance(measure, dGamma):
        return Core(set([(integrand, measure)]))
    else:
        raise RuntimeError('Illegal measure')


class Core(object):
    def __init__(self, set_c):
        self.cores = set_c

    def __add__(self, other):
        return Core(self.cores.union(other.cores))

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
