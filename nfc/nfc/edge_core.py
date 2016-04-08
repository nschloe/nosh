# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .discretize_edge_integral import DiscretizeEdgeIntegral
from .code_generator_eigen import CodeGeneratorEigen
from .helpers import extract_c_expression, templates_dir, is_affine_linear


def get_code_edge_core(namespace, class_name, core):
    '''Get code discretizer from raw core object.
    '''
    x0 = sympy.MatrixSymbol('x0', 3, 1)
    x1 = sympy.MatrixSymbol('x1', 3, 1)
    edge_length = sympy.Symbol('edge_length')
    edge_covolume = sympy.Symbol('edge_covolume')
    all_symbols = set([x0, x1, edge_length, edge_covolume])

    specs = inspect.getargspec(method)
    assert(len(specs.args) == len(all_symbols) + 1)

    edge_coeff, edge_affine = core.eval(x0, x1, edge_length, edge_covolume)

    return CodeEdgeCore(
            namespace, class_name,
            edge_coeff, edge_affine
            )


def get_edge_core_code_from_integral(
        namespace, class_name, u, integrand, measure
        ):
    '''Get code discretizer from expression.
    '''
    edge_coeff, edge_affine = _get_expressions_from_integrand(u, integrand)

    return _get_edge_core_code(
            namespace, class_name,
            edge_coeff, edge_affine,
            measure.subdomains
            )


def _get_expressions_from_integrand(u, integrand):
    x = sympy.MatrixSymbol('x', 3, 1)

    discretizer = DiscretizeEdgeIntegral()
    expr = discretizer.generate(integrand(x), u, x)

    if not is_affine_linear(expr, [discretizer.u0, discretizer.u1]):
        raise RuntimeError((
            'The given function\n'
            '    f(x) = %s\n'
            'does not seem to be affine linear in u.')
            % function(x)
            )

    # Get the coefficients of u0, u1.
    coeff00 = sympy.diff(expr, discretizer.u0)
    coeff01 = sympy.diff(expr, discretizer.u1)

    # Now construct the coefficients for the other way around.
    coeff10 = coeff01.subs([
        (discretizer.u0, discretizer.u1),
        (discretizer.u1, discretizer.u0),
        (nfl.n, nfl.neg_n)
        ])
    coeff11 = coeff00.subs([
        (discretizer.u0, discretizer.u1),
        (discretizer.u1, discretizer.u0),
        (nfl.n, nfl.neg_n)
        ])

    affine = expr.subs([(discretizer.u0, 0), (discretizer.u1, 0)])

    return (
        [[coeff00, coeff01], [coeff10, coeff11]],
        [affine, affine]
        )


def _get_edge_core_code(
        namespace, class_name,
        edge_coeff, edge_affine,
        subdomains
        ):
    assert(len(edge_coeff) == 2)
    assert(len(edge_coeff[0]) == 2)
    assert(len(edge_coeff[1]) == 2)
    assert(len(edge_affine) == 2)

    # edge
    edge_body, used_expressions = \
        _get_code_elements_edge(edge_coeff, edge_affine)

    # now take care of the template substitution
    members_init = []
    members_declare = []

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

    members_init.append(
            'nosh::edge_core(%s)' % parent_init
            )

    # init and declare expressions
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
    with open(os.path.join(templates_dir, 'edge_core.tpl'), 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': class_name,
            'edge00': _expr_to_code(edge_coeff[0][0]),
            'edge01': _expr_to_code(edge_coeff[0][1]),
            'edge10': _expr_to_code(edge_coeff[1][0]),
            'edge11': _expr_to_code(edge_coeff[1][1]),
            'edge_affine0': _expr_to_code(-edge_affine[0]),
            'edge_affine1': _expr_to_code(-edge_affine[1]),
            'edge_body': '\n'.join(edge_body),
            'members_init': members_init_code,
            'members_declare': '\n'.join(members_declare)
            })
    return code, dependencies


def _expr_to_code(expr):
    gen = CodeGeneratorEigen()
    code, _ = gen.generate(expr)
    return code


def _get_code_elements_edge(edge_coeff, edge_affine):
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
