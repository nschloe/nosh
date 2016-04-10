# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
import sympy

from .matrix_core_edge import get_matrix_core_edge_code_from_integral
from .matrix_core_vertex import get_matrix_core_vertex_code_from_integral
from .matrix_core_boundary import get_matrix_core_boundary_code_from_integral
from .helpers import templates_dir


def get_code_fvm_matrix(namespace, class_name, obj):
    if getattr(obj, 'dirichlet_boundary_conditions', None):
        dbcs = set(obj.dirichlet_boundary_conditions)
        for dbc in dbcs:
            assert(isinstance(dbc, nfl.DirichletBC))
    else:
        dbcs = set()

    dependencies = dbcs.copy()

    u = sympy.Function('u')
    res = obj.apply(u)

    matrix_core_edge_names = set()
    matrix_core_vertex_names = set()
    matrix_core_boundary_names = set()
    assert(isinstance(res, nfl.Core))
    code = ''
    for core in res.cores:
        integrand, measure, subdomains = core
        if isinstance(measure, nfl.dS):
            core_class_name = 'matrix_core_edge%d' % len(matrix_core_edge_names)
            core_code, deps = get_matrix_core_edge_code_from_integral(
                    namespace, core_class_name, u, integrand, subdomains
                    )
            matrix_core_edge_names.add(core_class_name)
        elif isinstance(measure, nfl.dV):
            core_class_name = 'matrix_core_vertex%d' % len(matrix_core_vertex_names)
            core_code, deps = get_matrix_core_vertex_code_from_integral(
                    namespace, core_class_name, u, integrand, subdomains
                    )
            matrix_core_vertex_names.add(core_class_name)
        elif isinstance(measure, nfl.dGamma):
            core_class_name = 'matrix_core_boundary%d' % len(matrix_core_boundary_names)
            core_code, deps = get_matrix_core_boundary_code_from_integral(
                    namespace, core_class_name, u, integrand, subdomains
                    )
            matrix_core_boundary_names.add(core_class_name)
        else:
            raise RuntimeError('Illegal measure type \'%s\'.' % measure)

        # since this object contains the cores, its dependencies contain the
        # dependencies of the cores as well
        dependencies.update(deps)

        code += '\n' + core_code

    # append the code of the linear problem itself
    code += '\n' + _get_code_linear_problem(
            class_name,
            matrix_core_edge_names,
            matrix_core_vertex_names,
            matrix_core_boundary_names,
            dbcs
            )

    return code, dependencies


def _get_code_linear_problem(
        name,
        matrix_core_edge_names,
        matrix_core_vertex_names,
        matrix_core_boundary_names,
        dbcs
        ):
    constructor_args = [
        'const std::shared_ptr<const nosh::mesh> & _mesh'
        ]
    init_matrix_core_edges = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in matrix_core_edge_names]
                )
            )
    init_matrix_core_vertexs = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in matrix_core_vertex_names]
                )
            )
    init_matrix_core_boundarys = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in matrix_core_boundary_names]
                )
            )

    # handle the boundary conditions
    joined = ', '.join(
            'std::make_shared<%s>()' %
            type(bc).__name__.lower() for bc in dbcs
            )
    init_dbcs = '{%s}' % joined
    # boundary conditions handling done

    members_init = [
      'nosh::fvm_matrix(\n_mesh,\n %s,\n %s,\n %s,\n %s\n)' %
      (init_matrix_core_edges, init_matrix_core_vertexs, init_matrix_core_boundarys, init_dbcs)
      ]
    members_declare = []

    templ = os.path.join(templates_dir, 'fvm_matrix.tpl')
    with open(templ, 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': name.lower(),
            'constructor_args': ',\n'.join(constructor_args),
            'members_init': ',\n'.join(members_init),
            'members_declare': '\n'.join(members_declare)
            })

    return code
