# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
import sympy

# from .operator_core_boundary \
#     import get_operator_core_boundary_code_from_integral
from .operator_core_dirichlet import get_code_dirichlet
# from .operator_core_edge import get_operator_core_edge_code_from_integral
from .operator_core_vertex import get_operator_core_vertex_code_from_integral
from .helpers import get_uuid, templates_dir


def get_code_fvm_operator(namespace, class_name, obj):

    # code, dependencies, operator_core_names, fvm_operator_names = \
    #        handle_dependencies(namespace, obj)

    code, dependencies, operator_core_names = \
            handle_dependencies(namespace, obj)

    # append the code of the linear problem itself
    code += '\n' + get_code_operator(
            'fvm_operator.tpl',
            class_name.lower(),
            'nosh::fvm_operator',
            operator_core_names['edge'],
            operator_core_names['vertex'],
            operator_core_names['boundary'],
            operator_core_names['dirichlet'],
            operator_core_names['operator']
            )

    return code, dependencies


def handle_dependencies(namespace, obj):
    dependencies = set()

    u = sympy.Function('u')
    res = obj.apply(u)

    assert(isinstance(res, nfl.CoreList))

    integrals_code, integrals_deps, core_names = \
        handle_integrals(namespace, u, res.integrals)

    dirichlet_code, dirichlet_deps, dcore_names = \
        handle_dirichlets(namespace, obj.dirichlet)

    core_names['dirichlet'] = dcore_names

    operators_deps = set(res.fvm_matrices)
    core_names['operator'] = \
        set([op.__class__.__name__.lower() for op in res.fvm_matrices])

    code = '\n'.join([integrals_code, dirichlet_code])
    dependencies = set().union(
            integrals_deps,
            dirichlet_deps,
            operators_deps
            )

    return code, dependencies, core_names


def handle_integrals(namespace, u, integrals):
    operator_core_names = {
            'edge': set(),
            'vertex': set(),
            'boundary': set(),
            'dirichlet': [],
            }

    core_codes = []
    dependencies = set()

    for integral in integrals:
        if isinstance(integral.measure, nfl.ControlVolumeSurface):
            core_class_name = 'operator_core_edge_%s' % get_uuid()
            core_code, deps = get_operator_core_edge_code_from_integral(
                    namespace, core_class_name, u,
                    integral.integrand, integral.subdomains
                    )
            operator_core_names['edge'].add(core_class_name)
        elif isinstance(integral.measure, nfl.ControlVolume):
            core_class_name = 'operator_core_vertex_%s' % get_uuid()
            core_code, deps = get_operator_core_vertex_code_from_integral(
                    namespace, core_class_name, u,
                    integral.integrand, integral.subdomains
                    )
            operator_core_names['vertex'].add(core_class_name)
        elif isinstance(integral.measure, nfl.BoundarySurface):
            core_class_name = 'operator_core_boundary_%s' % get_uuid()
            core_code, deps = get_operator_core_boundary_code_from_integral(
                    namespace, core_class_name, u,
                    integral.integrand, integral.subdomains
                    )
            operator_core_names['boundary'].add(core_class_name)
        else:
            raise RuntimeError('Illegal measure type \'%s\'.' % measure)

        # since this object contains the cores, its dependencies contain the
        # dependencies of the cores as well
        dependencies.update(deps)
        core_codes.append(core_code)

    return '\n'.join(core_codes), dependencies, operator_core_names


def handle_dirichlets(namespace, dirichlets):
    names = []
    for k, dirichlet in enumerate(dirichlets):
        f, subdomains = dirichlet

        if not isinstance(subdomains, list):
            try:
                subdomains = list(subdomains)
            except TypeError:  # TypeError: 'D1' object is not iterable
                subdomains = [subdomains]

        core_class_name = 'operator_core_dirichlet_%s' % get_uuid()
        core_code, dependencies = \
            get_code_dirichlet(core_class_name, f, subdomains)
        names.append(core_class_name)

    return core_code, dependencies, names


def get_code_operator(
        template_filename,
        class_name,
        base_class_name,
        core_edge_names,
        core_vertex_names,
        core_boundary_names,
        core_dirichlet_names,
        core_operator_names
        ):
    constructor_args = [
        'const std::shared_ptr<const nosh::mesh> & _mesh'
        ]
    init_core_edge = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n for n in core_edge_names]
                )
            )
    init_core_vertex = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n for n in core_vertex_names]
                )
            )
    init_core_boundary = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n for n in core_boundary_names]
                )
            )
    init_core_dirichlet = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n for n in core_dirichlet_names]
                )
            )
    init_operator = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>(_mesh)' % n for n in core_operator_names]
                )
            )

    members_init = [
      '%s(\n_mesh,\n %s,\n %s,\n %s,\n %s,\n %s\n)' %
      (base_class_name,
       init_core_edge,
       init_core_vertex,
       init_core_boundary,
       init_core_dirichlet,
       init_operator
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
