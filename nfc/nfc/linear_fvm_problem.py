# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
import sympy

from .edge_core import get_edge_core_code_from_integrand
from .vertex_core import get_vertex_core_code_from_integrand
from .boundary_core import get_boundary_core_code_from_integrand
from .helpers import templates_dir


def get_code_linear_fvm_problem(namespace, class_name, obj):
    if getattr(obj, 'dirichlet_boundary_conditions', None):
        dbcs = set(obj.dirichlet_boundary_conditions)
        for dbc in dbcs:
            assert(isinstance(dbc, nfl.DirichletBC))
    else:
        dbcs = set()

    dependencies = dbcs.copy()

    u = sympy.Function('u')
    res = obj.eval(u)

    edge_core_names = set()
    vertex_core_names = set()
    boundary_core_names = set()
    assert(isinstance(res, nfl.Core))
    code = ''
    for core in res.cores:
        integrand, measure = core
        if measure == 'dS':
            core_class_name = 'edge_core%d' % len(edge_core_names)
            core_code, deps = get_edge_core_code_from_integrand(
                    namespace, core_class_name, u, integrand
                    )
            edge_core_names.add(core_class_name)
        elif measure == 'dV':
            core_class_name = 'vertex_core%d' % len(vertex_core_names)
            core_code, deps = get_vertex_core_code_from_integrand(
                    namespace, core_class_name, u, integrand
                    )
            vertex_core_names.add(core_class_name)
        elif measure == 'dGamma':
            core_class_name = 'boundary_core%d' % len(boundary_core_names)
            core_code, deps = get_boundary_core_code_from_integrand(
                    namespace, core_class_name, u, integrand
                    )
            boundary_core_names.add(core_class_name)
        else:
            raise RuntimeError('Illegal measure type \'%s\'.' % measure)

        # since this object contains the cores, its dependencies contain the
        # dependencies of the cores as well
        dependencies.update(deps)

        code += '\n' + core_code

    # append the code of the linear problem itself
    code += '\n' + _get_code_linear_problem(
            class_name,
            edge_core_names,
            vertex_core_names,
            boundary_core_names,
            dbcs
            )

    return code, dependencies


def _get_code_linear_problem(
        name,
        edge_core_names,
        vertex_core_names,
        boundary_core_names,
        dbcs
        ):
    constructor_args = [
        'const std::shared_ptr<const nosh::mesh> & _mesh'
        ]
    init_edge_cores = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in edge_core_names]
                )
            )
    init_vertex_cores = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in vertex_core_names]
                )
            )
    init_boundary_cores = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in boundary_core_names]
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
      'nosh::linear_problem(\n_mesh,\n %s,\n %s,\n %s,\n %s\n)' %
      (init_edge_cores, init_vertex_cores, init_boundary_cores, init_dbcs)
      ]
    members_declare = []

    templ = os.path.join(templates_dir, 'linear_fvm_problem.tpl')
    with open(templ, 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': name.lower(),
            'constructor_args': ',\n'.join(constructor_args),
            'members_init': ',\n'.join(members_init),
            'members_declare': '\n'.join(members_declare)
            })

    return code
