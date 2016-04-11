# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .code_generator_eigen import CodeGeneratorEigen
from .helpers import \
        compare_variables, \
        extract_c_expression, \
        templates_dir, \
        is_affine_linear


def get_operator_core_vertex_code_from_integral(
        namespace, class_name,
        u, integrand, subdomains
        ):
    '''Get code generator from expression.
    '''
    expr = _get_expression_from_integral(u, integrand)

    return _get_code_operator_core_vertex(
            namespace, class_name,
            expr,
            subdomains
            )


def _get_expression_from_integral(u, function):
    # Numerically integrate function over a control volume.
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
    # Make sure that f is affine linear in u0; we're building a operator
    # here.
    control_volume = sympy.Symbol('control_volume')
    return control_volume * fu0


def _get_code_operator_core_vertex(
        namespace, class_name,
        expr,
        subdomains
        ):
    arguments = set([
        sympy.MatrixSymbol('x', 3, 1),
        sympy.Symbol('control_volume'),
        sympy.Symbol('u0')
        ])
    unused_arguments, used_expressions = compare_variables(arguments, [expr])

    # now take care of the template substitution
    members_init = []
    members_declare = []

    # init and declare all expressions
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
    members_init.append('nosh::operator_core_vertex(%s)' % parent_init)

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
    filename = os.path.join(templates_dir, 'operator_core_vertex.tpl')
    with open(filename, 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': class_name,
            'return_value': extract_c_expression(expr),
            'vertex_body': '\n'.join(
                ('(void) %s;' % name) for name in unused_arguments
                ),
            'members_init': members_init_code,
            'members_declare': '\n'.join(members_declare)
            })
    return code, dependencies
