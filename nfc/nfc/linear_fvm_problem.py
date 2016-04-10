# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
import sympy

from .matrix_core_boundary import get_matrix_core_boundary_code_from_integral
from .matrix_core_dirichlet import get_code_dirichlet
from .matrix_core_edge import get_matrix_core_edge_code_from_integral
from .matrix_core_vertex import get_matrix_core_vertex_code_from_integral
from .helpers import templates_dir


def get_code_linear_fvm_problem(namespace, class_name, obj):

    dependencies = set()

    u = sympy.Function('u')
    res = obj.eval(u)

    matrix_core_edge_names = set()
    matrix_core_vertex_names = set()
    matrix_core_boundary_names = set()
    assert(isinstance(res, nfl.Core))
    code = ''
    for core in res.cores:
        integrand, measure, subdomains = core
        if isinstance(measure, nfl.dS):
            core_class_name = \
              'matrix_core_edge%d' % len(matrix_core_edge_names)
            core_code, deps = get_matrix_core_edge_code_from_integral(
                    namespace, core_class_name, u, integrand, subdomains
                    )
            matrix_core_edge_names.add(core_class_name)
        elif isinstance(measure, nfl.dV):
            core_class_name = \
                'matrix_core_vertex%d' % len(matrix_core_vertex_names)
            core_code, deps = get_matrix_core_vertex_code_from_integral(
                    namespace, core_class_name, u, integrand, subdomains
                    )
            matrix_core_vertex_names.add(core_class_name)
        elif isinstance(measure, nfl.dGamma):
            core_class_name = \
                'matrix_core_boundary%d' % len(matrix_core_boundary_names)
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

    matrix_core_dirichlet_names = []
    for k, dirichlet in enumerate(obj.dirichlet):
        f, subdomains = dirichlet

        if not isinstance(subdomains, list):
            try:
                subdomains = list(subdomains)
            except TypeError:  # TypeError: 'D1' object is not iterable
                subdomains = [subdomains]

        core_class_name = 'matrix_core_dirichlet%d' % k
        core_code, deps = get_code_dirichlet(core_class_name, f, subdomains)
        matrix_core_dirichlet_names.append(core_class_name)
        dependencies.update(deps)
        code += '\n' + core_code

    # append the code of the linear problem itself
    code += '\n' + _get_code_linear_problem(
            class_name,
            matrix_core_edge_names,
            matrix_core_vertex_names,
            matrix_core_boundary_names,
            matrix_core_dirichlet_names
            )

    return code, dependencies


def _get_code_linear_problem(
        name,
        matrix_core_edge_names,
        matrix_core_vertex_names,
        matrix_core_boundary_names,
        matrix_core_dirichlet_names
        ):
    constructor_args = [
        'const std::shared_ptr<const nosh::mesh> & _mesh'
        ]
    init_matrix_core_edge = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in matrix_core_edge_names]
                )
            )
    init_matrix_core_vertex = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in matrix_core_vertex_names]
                )
            )
    init_matrix_core_boundary = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in matrix_core_boundary_names]
                )
            )
    init_matrix_core_dirichlet = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in matrix_core_dirichlet_names]
                )
            )

    members_init = [
      'nosh::linear_problem(\n_mesh,\n %s,\n %s,\n %s,\n %s\n)' %
      (init_matrix_core_edge, init_matrix_core_vertex,
          init_matrix_core_boundary, init_matrix_core_dirichlet)
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
