# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .code_generator_eigen import CodeGeneratorEigen
from .helpers import extract_c_expression, templates_dir, is_affine_linear


def get_vertex_core_code(namespace, class_name, core):
    '''Get code generator from raw core object.
    '''
    # handle the vertex contributions
    x = sympy.MatrixSymbol('x')
    vol = sympy.Symbol('control_volume')
    all_symbols = set([x, vol])

    specs = inspect.getargspec(method)
    assert(len(specs.args) == len(all_symbols) + 1)

    vertex_coeff, vertex_affine = method(x, vol)

    return _get_code_vertex_core(
            namespace, class_name,
            vertex_coeff, vertex_affine
            )


def get_vertex_core_code_from_integral(
        namespace, class_name,
        u, integrand, subdomains
        ):
    '''Get code generator from expression.
    '''
    vertex_coeff, vertex_affine = _get_expressions_from_integral(u, integrand)

    return _get_code_vertex_core(
            namespace, class_name,
            vertex_coeff, vertex_affine,
            subdomains
            )


def _get_expressions_from_integral(u, function):
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
    assert(is_affine_linear(fu0, [u0]))
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


def _get_code_vertex_core(
        namespace, class_name,
        vertex_coeff, vertex_affine,
        subdomains
        ):
    # vertex
    vertex_arguments = set([
        sympy.MatrixSymbol('x', 3, 1),
        sympy.Symbol('control_volume')
        ])
    vertex_unused_arguments, vertex_used_expressions = \
        _scan_code(vertex_arguments, [vertex_coeff, vertex_affine])

    # now take care of the template substitution
    members_init = []
    members_declare = []

    # init and declare all expressions
    used_expressions = vertex_used_expressions

    dependencies = set()
    dependencies.update(used_expressions)

    # init parent object
    subdomain_ids = set([
        sd.__class__.__name__.lower() for sd in subdomains
        ])
    if len(subdomain_ids) == 0:
        # If nothing is specified, use the entire boundary
        subdomain_ids.add('everywhere')
    parent_init = '{%s}' % ', '.join(['"%s"' % s for s in subdomain_ids])
    members_init.append('nosh::vertex_core(%s)' % parent_init)

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
    with open(os.path.join(templates_dir, 'vertex_core.tpl'), 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': class_name,
            'vertex_contrib': extract_c_expression(vertex_coeff),
            'vertex_affine': extract_c_expression(-vertex_affine),
            'vertex_body': '\n'.join(
                ('(void) %s;' % name) for name in vertex_unused_arguments
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
