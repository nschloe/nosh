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


def get_code_fvm_matrix(namespace, class_name, obj):

    code, dependencies, matrix_core_names = \
            handle_core_dependencies(namespace, obj)

    # append the code of the linear problem itself
    code += '\n' + get_code_linear_problem(
            'fvm_matrix.tpl',
            class_name.lower(),
            'nosh::fvm_matrix',
            matrix_core_names['edge'],
            matrix_core_names['vertex'],
            matrix_core_names['boundary'],
            matrix_core_names['dirichlet'],
            )

    return code, dependencies


def handle_core_dependencies(namespace, obj):
    dependencies = set()

    u = sympy.Function('u')
    res = obj.apply(u)

    assert(isinstance(res, nfl.CoreList))

    matrix_core_names = {
            'edge': set(),
            'vertex': set(),
            'boundary': set(),
            'dirichlet': []
            }

    code = ''
    for integral in res.integrals:
        if isinstance(integral.measure, nfl.dS):
            core_class_name = \
              'matrix_core_edge%d' % len(matrix_core_names['edge'])
            core_code, deps = get_matrix_core_edge_code_from_integral(
                    namespace, core_class_name, u,
                    integral.integrand, integral.subdomains
                    )
            matrix_core_names['edge'].add(core_class_name)
        elif isinstance(integral.measure, nfl.dV):
            core_class_name = \
                'matrix_core_vertex%d' % len(matrix_core_names['vertex'])
            core_code, deps = get_matrix_core_vertex_code_from_integral(
                    namespace, core_class_name, u,
                    integral.integrand, integral.subdomains
                    )
            matrix_core_names['vertex'].add(core_class_name)
        elif isinstance(integral.measure, nfl.dGamma):
            core_class_name = \
                'matrix_core_boundary%d' % len(matrix_core_names['boundary'])
            core_code, deps = get_matrix_core_boundary_code_from_integral(
                    namespace, core_class_name, u,
                    integral.integrand, integral.subdomains
                    )
            matrix_core_names['boundary'].add(core_class_name)
        else:
            raise RuntimeError('Illegal measure type \'%s\'.' % measure)

        # since this object contains the cores, its dependencies contain the
        # dependencies of the cores as well
        dependencies.update(deps)
        code += '\n' + core_code

    for k, dirichlet in enumerate(obj.dirichlet):
        f, subdomains = dirichlet

        if not isinstance(subdomains, list):
            try:
                subdomains = list(subdomains)
            except TypeError:  # TypeError: 'D1' object is not iterable
                subdomains = [subdomains]

        core_class_name = 'matrix_core_dirichlet%d' % k
        core_code, deps = get_code_dirichlet(core_class_name, f, subdomains)
        matrix_core_names['dirichlet'].append(core_class_name)
        dependencies.update(deps)
        code += '\n' + core_code

    return code, dependencies, matrix_core_names


def get_code_linear_problem(
        template_filename,
        class_name,
        base_class_name,
        matrix_core_edge_names,
        matrix_core_vertex_names,
        matrix_core_boundary_names,
        matrix_core_dirichlet_names,
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
      '%s(\n_mesh,\n %s,\n %s,\n %s,\n %s\n)' %
      (base_class_name,
       init_matrix_core_edge,
       init_matrix_core_vertex,
       init_matrix_core_boundary,
       init_matrix_core_dirichlet
       )
      ]
    members_declare = []

    templ = os.path.join(templates_dir, template_filename)
    with open(templ, 'r') as f:
        src = Template(f.read())
        code = src.substitute({
            'name': class_name,
            'constructor_args': ',\n'.join(constructor_args),
            'members_init': ',\n'.join(members_init),
            'members_declare': '\n'.join(members_declare)
            })

    return code
