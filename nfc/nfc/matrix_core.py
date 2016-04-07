# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .discretize_edge_integral import DiscretizeEdgeIntegral
from .code_generator_eigen import CodeGeneratorEigen
from .helpers import extract_c_expression, templates_dir


def _is_affine_linear(expr, vars):
    for var in vars:
        if not sympy.Eq(sympy.diff(expr, var, var), 0):
            return False
    return True


# We still need this for pure matrices
# def _is_linear(expr, vars):
#     if not _is_affine_linear(expr, vars):
#         return False
#     # Check that expr is not affine.
#     if isinstance(expr, int) or isinstance(expr, float):
#         return expr == 0
#     else:
#         return expr.subs([(var, 0) for var in vars]) == 0


def _get_expr_edge_contrib(self, method):
    x0 = sympy.MatrixSymbol('x0', 3, 1)
    x1 = sympy.MatrixSymbol('x1', 3, 1)
    edge_length = sympy.Symbol('edge_length')
    edge_covolume = sympy.Symbol('edge_covolume')
    all_symbols = set([x0, x1, edge_length, edge_covolume])

    specs = inspect.getargspec(method)
    assert(len(specs.args) == len(all_symbols) + 1)

    result = method(x0, x1, edge_length, edge_covolume)

    assert 'coeff' in result
    assert(len(result['coeff']) == 2)
    assert(len(result['coeff'][0]) == 2)
    assert(len(result['coeff'][1]) == 2)
    assert 'affine' in result
    assert(len(result['affine']) == 2)

    return result['coeff'], result['affine']


def _get_expr_vertex_contrib(self, method):
    # handle the vertex contributions
    x = sympy.MatrixSymbol('x')
    vol = sympy.Symbol('control_volume')
    all_symbols = set([x, vol])

    specs = inspect.getargspec(method)
    assert(len(specs.args) == len(all_symbols) + 1)

    result = method(x, vol)

    return result['coeff'], result['affine']


def get_core_code_generator(namespace, name, core):
    '''Get code generator from raw core object.
    '''
    # handle the edge contributions
    if callable(getattr(core, 'edge_contrib', None)):
        edge_coeff, edge_affine = \
                _get_expr_edge_contrib(core.edge_contrib)

    if callable(getattr(core, 'vertex_contrib', None)):
        vertex_coeff, vertex_affine = \
                _get_expr_vertex_contrib(core.vertex_contrib)

    if callable(getattr(core, 'domain_boundary_contrib', None)):
        db_coeff, db_affine = \
                _get_expr_db_contrib(core.vertex_contrib)
    return CodeMatrixCore(
            namespace, name,
            edge_coeff, edge_affine,
            vertex_coeff, vertex_affine,
            db_coeff, db_affine
            )


def get_core_code_from_expression(namespace, name, expr):
    '''Get code generator from expression.
    '''
    u = sympy.Function('u')
    res = expr(u)
    assert(isinstance(res, nfl.Core))

    edge_coeff, edge_affine = _get_expr_edge(u, res.edge)
    vertex_coeff, vertex_affine = _get_expr_vertex(u, res.vertex)
    db_coeff, db_affine = _get_expr_db(res.domain_boundary)

    return CodeMatrixCore(
            namespace, name,
            edge_coeff, edge_affine,
            vertex_coeff, vertex_affine,
            db_coeff, db_affine
            )


def _get_expr_edge(u, function):
    x = sympy.MatrixSymbol('x', 3, 1)

    generator = DiscretizeEdgeIntegral()
    expr, _ = generator.generate(function(x), u, x)

    # Make sure the expression is linear in u0, u1.
    if not _is_affine_linear(expr, [generator.u0, generator.u1]):
        raise RuntimeError((
            'The given function\n'
            '    f(x) = %s\n'
            'does not seem to be affine linear in u.')
            % function(x)
            )

    # Get the coefficients of u0, u1.
    coeff00 = sympy.diff(expr, generator.u0)
    coeff01 = sympy.diff(expr, generator.u1)

    # Now construct the coefficients for the other way around.
    coeff10 = coeff01.subs([
        (generator.u0, generator.u1),
        (generator.u1, generator.u0),
        (nfl.n, nfl.neg_n)
        ])
    coeff11 = coeff00.subs([
        (generator.u0, generator.u1),
        (generator.u1, generator.u0),
        (nfl.n, nfl.neg_n)
        ])

    affine = expr.subs([(generator.u0, 0), (generator.u1, 0)])

    return (
        [[coeff00, coeff01], [coeff10, coeff11]],
        [affine, affine]
        )


def _get_expr_vertex(u, function):
    # Numerically integrate function over a control volume.
    # TODO find out if the code can be evaluated at the volume "midpoint",
    # then evaluate it there, and multiply by the volume.
    x = sympy.MatrixSymbol('x', 3, 1)
    # Evaluate the function for u at x.
    fx = function(x)
    # Replace all occurences of u(x) by u0 (the value at the control volume
    # center) and multiply by the control volume)
    u0 = sympy.Symbol('u0')
    try:
        fu0 = fx.subs(u(x), u0)
    except AttributeError:
        # 'int' object has no attribute 'subs'
        fu0 = fx
    # Make sure that f is affine linear in u0; we're building a matrix
    # here.
    assert(_is_affine_linear(fu0, [u0]))
    control_volume = sympy.Symbol('control_volume')
    # Get coefficient of u0
    coeff = sympy.diff(fu0, u0)
    coeff = control_volume * coeff
    # Get affine part
    if isinstance(fu0, float):
        affine = control_volume * fu0
    else:
        affine = control_volume * fu0.subs(u0, 0)
    return coeff, affine


def _get_expr_db(function):
    # Numerically integrate function over a part of the domain boundary,
    # namely the part surrounding a vertex.
    x = sympy.MatrixSymbol('x', 3, 1)
    # Evaluate the function for u at x.
    fx = function(x)
    # Replace all occurences of u(x) by u0 (the value at the center) and
    # multiply by the surface area)
    u0 = sympy.Symbol('u0')
    try:
        fu0 = fx.subs(u(x), u0)
    except AttributeError:
        # 'int' object has no attribute 'subs'
        fu0 = fx
    # Make sure that f is affine linear in u0; we're building a matrix
    # here.
    assert(_is_affine_linear(fu0, [u0]))
    control_volume = sympy.Symbol('surface_area')
    # Get coefficient of u0
    coeff = sympy.diff(fu0, u0)
    coeff = control_volume * coeff
    # Get affine part
    if isinstance(fu0, float) or isinstance(fu0, int):
        affine = control_volume * fu0
    else:
        affine = control_volume * fu0.subs(u0, 0)
    return (coeff, affine)


class CodeMatrixCore(object):
    def __init__(
            self, namespace, name,
            edge_coeff, edge_affine,
            vertex_coeff, vertex_affine,
            db_coeff, db_affine
            ):
        # edge
        edge_body, edge_used_expressions = \
            self.get_code_elements_edge(edge_coeff, edge_affine)

        # vertex
        vertex_arguments = set([
            sympy.MatrixSymbol('x', 3, 1),
            sympy.Symbol('control_volume')
            ])
        vertex_unused_arguments, vertex_used_expressions = \
            self.scan_code(vertex_arguments, [vertex_coeff, vertex_affine])

        # domain boundary
        db_arguments = set([
            sympy.MatrixSymbol('x', 3, 1),
            sympy.Symbol('surface_area')
            ])
        db_unused_arguments, db_used_expressions = \
            self.scan_code(db_arguments, [db_coeff, db_affine])

        # now take care of the template substitution
        members_init = []
        members_declare = []

        # init and declare all expressions
        used_expressions = set().union(
                edge_used_expressions,
                vertex_used_expressions,
                db_used_expressions
                )

        self.dependencies = set()
        self.dependencies.update(used_expressions)

        # init and declare expressions in the C++ code
        for expr in used_expressions:
            members_init.append('%s(%s::%s())' % (expr, self.namespace, expr))
            members_declare.append(
                    'const %s::%s %s;' % (self.namespace, expr, expr)
                    )

        if members_init:
            members_init_code = ':\n' + ',\n'.join(members_init)
        else:
            members_init_code = ''

        # template substitution
        with open(os.path.join(templates_dir, 'matrix_core.tpl'), 'r') as f:
            src = Template(f.read())
            self.code = src.substitute({
                'name': name.lower() + '_core',
                'edge00': self.expr_to_code(edge_coeff[0][0]),
                'edge01': self.expr_to_code(edge_coeff[0][1]),
                'edge10': self.expr_to_code(edge_coeff[1][0]),
                'edge11': self.expr_to_code(edge_coeff[1][1]),
                'edge_affine0': self.expr_to_code(-edge_affine[0]),
                'edge_affine1': self.expr_to_code(-edge_affine[1]),
                'edge_body': '\n'.join(edge_body),
                'vertex_contrib': extract_c_expression(vertex_coeff),
                'vertex_affine': extract_c_expression(-vertex_affine),
                'vertex_body': '\n'.join(
                    ('(void) %s;' % name) for name in vertex_unused_arguments
                    ),
                'db_coeff': extract_c_expression(db_coeff),
                'db_affine': extract_c_expression(-db_affine),
                'db_body': '\n'.join(
                    ('(void) %s;' % name) for name in db_unused_arguments
                    ),
                'members_init': members_init_code,
                'members_declare': '\n'.join(members_declare)
                })
        return

    def get_dependencies(self):
        return self.dependencies

    def get_code(self):
        return self.code

    def expr_to_code(self, expr):
        gen = CodeGeneratorEigen()
        code, _ = gen.generate(expr)
        return code

    def get_code_elements_edge(self, edge_coeff, edge_affine):
        # Check out used an unused symbols for the edge code.
        # Arguments from the template.
        arguments = set([
            sympy.Symbol('x0'),
            sympy.Symbol('x1'),
            sympy.Symbol('edge_length'),
            sympy.Symbol('edge_covolume')
            ])
        expressions = [
            edge_coeff[0][0], edge_coeff[0][1],
            edge_coeff[1][0], edge_coeff[1][1],
            edge_affine[0], edge_affine[1]
            ]

        used_symbols = set([])
        used_expressions = set([])
        body = []
        for expr in expressions:
            try:
                used_symbols.update(expr.free_symbols)
            except AttributeError:
                pass

            used_expressions.update(set([
                    type(atom) for atom in expr.atoms(nfl.Expression)
                    ]))

        unused_arguments = arguments - used_symbols
        undefined_symbols = used_symbols - arguments

        # special treatment for n
        if nfl.n in undefined_symbols or \
           nfl.neg_n in undefined_symbols:
            body.append('const auto n = (x1 - x0) / edge_length;')
            unused_arguments -= set([
                sympy.Symbol('x0'),
                sympy.Symbol('x1')
                ])
            undefined_symbols.remove(nfl.n)
            undefined_symbols.remove(nfl.neg_n)

        assert(len(undefined_symbols) == 0)

        for name in unused_arguments:
            body.append('(void) %s;' % name)

        return body, used_expressions

    def scan_code(self, arguments, expressions):
        used_symbols = set([])
        used_expressions = set([])

        for expr in expressions:
            try:
                used_symbols.update(expr.free_symbols)
            except AttributeError:
                pass

            used_expressions.update(set([
                    type(atom) for atom in expr.atoms(nfl.Expression)
                    ]))

        unused_arguments = arguments - used_symbols
        undefined_symbols = used_symbols - arguments
        assert(len(undefined_symbols) == 0)

        return unused_arguments, used_expressions
