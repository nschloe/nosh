# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
import sympy
from .helpers import extract_c_expression, templates_dir


def get_code_dirichletbc(name, dbc):
    x = sympy.MatrixSymbol('x', 3, 1)

    dependencies = dbc.subdomains

    subdomain_ids = set([
        sd.__class__.__name__.lower() for sd in dbc.subdomains
        ])

    result = dbc.eval(x)
    try:
        is_x_used_eval = x in result.free_symbols
    except AttributeError:
        is_x_used_eval = False

    # template substitution
    with open(os.path.join(templates_dir, 'dirichlet_bc.tpl'), 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': name.lower(),
            'init': 'nosh::dirichlet_bc({%s})' % ', '.join([
                '"%s"' % s for s in subdomain_ids
                ]),
            'eval_return_value': extract_c_expression(result),
            'eval_void': '' if is_x_used_eval else '(void) x;\n',
            })

    return code, dependencies
