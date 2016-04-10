# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .discretize_edge_integral import DiscretizeEdgeIntegral
from .code_generator_eigen import CodeGeneratorEigen
from .helpers import extract_c_expression, templates_dir, is_affine_linear


def get_matrix_core_boundary_code_from_integral(
        namespace, class_name,
        u, integrand, subdomains
        ):
    '''Get code generator from expression.
    '''
    coeff, affine = _get_expressions_from_integral(integrand)

    return _get_code_matrix_core_boundary(
            namespace, class_name,
            coeff, affine,
            subdomains
            )


def _get_expressions_from_integral(function):
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
    assert(is_affine_linear(fu0, [u0]))
    control_volume = sympy.Symbol('surface_area')
    # Get coefficient of u0
    coeff = sympy.diff(fu0, u0)
    coeff = control_volume * coeff
    # Get affine part
    if isinstance(fu0, float) or isinstance(fu0, int):
        affine = control_volume * fu0
    else:
        affine = control_volume * fu0.subs(u0, 0)
    return coeff, affine


def _get_code_matrix_core_boundary(
        namespace, class_name,
        db_coeff, db_affine,
        subdomains
        ):
    # domain boundary
    db_arguments = set([
        sympy.MatrixSymbol('x', 3, 1),
        sympy.Symbol('surface_area')
        ])
    db_unused_arguments, db_used_expressions = \
        _scan_code(db_arguments, [db_coeff, db_affine])

    # now take care of the template substitution
    members_init = []
    members_declare = []

    # init parent object
    subdomain_ids = set([
        sd.__class__.__name__.lower() for sd in subdomains
        ])
    if len(subdomain_ids) == 0:
        # If nothing is specified, use the entire boundary
        subdomain_ids.add('everywhere')
    parent_init = '{%s}' % ', '.join(['"%s"' % s for s in subdomain_ids])
    members_init.append('nosh::matrix_core_vertex(%s)' % parent_init)

    # init and declare all expressions
    used_expressions = db_used_expressions

    dependencies = set()
    dependencies.update(used_expressions)

    # init and declare expressions in the C++ code
    for expr in used_expressions:
        members_init.append('%s(%s::%s())' % (expr, namespace, expr))
        members_declare.append(
                'const %s::%s %s;' % (namespace, expr, expr)
                )

    if members_init:
        members_init_code = ':\n' + ',\n'.join(members_init)
    else:
        members_init_code = ''

    # template substitution
    with open(os.path.join(templates_dir, 'matrix_core_boundary.tpl'), 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': class_name,
            'db_coeff': extract_c_expression(db_coeff),
            'db_affine': extract_c_expression(-db_affine),
            'db_body': '\n'.join(
                ('(void) %s;' % name) for name in db_unused_arguments
                ),
            'members_init': members_init_code,
            'members_declare': '\n'.join(members_declare)
            })
    return code, dependencies


def _expr_to_code(expr):
    gen = CodeGeneratorEigen()
    code, _ = gen.generate(expr)
    return code


def _scan_code(arguments, expressions):
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
