# -*- coding: utf-8 -*-
#
import inspect
import os
from string import Template
import sympy

import nfl

from .fvm_matrix import FvmMatrixCode, gather_core_dependencies
from .integral_boundary import IntegralBoundary
from .dirichlet import Dirichlet
from .integral_edge import IntegralEdge
from .integral_vertex import IntegralVertex
from .helpers import get_uuid, sanitize_identifier, templates_dir


class FvmOperatorCode(object):
    def __init__(self, namespace, cls):
        self.class_name = sanitize_identifier(cls.__name__)
        self.namespace = namespace
        self.dependencies = \
            gather_dependencies(namespace, cls, is_matrix=False)
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_object(self, dep_class_objects):
        # Go through the dependencies collect the cores.
        dirichlet_core_names = []
        vertex_core_names = []
        edge_core_names = []
        boundary_core_names = []
        fvm_matrix_names = []
        for dep in self.dependencies:
            if isinstance(dep, Dirichlet):
                dirichlet_core_names.append(dep.class_name)
            elif isinstance(dep, IntegralVertex):
                vertex_core_names.append(dep.class_name)
            elif isinstance(dep, IntegralEdge):
                edge_core_names.append(dep.class_name)
            elif isinstance(dep, IntegralBoundary):
                boundary_core_names.append(dep.class_name)
            elif isinstance(dep, FvmMatrixCode):
                fvm_matrix_names.append(dep.class_name)
            else:
                raise RuntimeError(
                    'Dependency \'%s\' not accounted for.' % dep.class_name
                    )

        code = get_code_linear_problem(
            'fvm_operator.tpl',
            self.class_name,
            'nosh::fvm_operator',
            edge_core_names,
            vertex_core_names,
            boundary_core_names,
            dirichlet_core_names,
            fvm_matrix_names
            )

        return {
            'code': code
            }


def gather_dependencies(namespace, cls, is_matrix):
    u = sympy.Function('u')
    u0 = sympy.Function('u0')

    if (len(inspect.getargspec(cls.apply).args) == 1):
        res = cls.apply(u)
    elif (len(inspect.getargspec(cls.apply).args) == 2):
        res = cls.apply(u, u0)
    else:
        raise ValueError('Only methods with one or two arguments allowed.')

    dependencies = gather_core_dependencies(namespace, expr, is_matrix)

    # Add dependencies on fvm_matrices
    for fvm_matrix in res.fvm_matrices:
        dependencies.add(
            FvmMatrixCode(namespace, fvm_matrix.__class__)
            )

    return dependencies


def get_code_fvm_operator(namespace, class_name, obj):
    code, dependencies, operator_core_names = \
            handle_core_dependencies(namespace, obj)

    # append the code of the linear problem itself
    code += '\n' + get_code_linear_problem(
            'fvm_operator.tpl',
            class_name.lower(),
            'nosh::fvm_operator',
            operator_core_names['edge'],
            operator_core_names['vertex'],
            operator_core_names['boundary'],
            operator_core_names['dirichlet'],
            )

    return code, dependencies


def get_code_linear_problem(
        template_filename,
        class_name,
        base_class_name,
        operator_core_edge_names,
        operator_core_vertex_names,
        operator_core_boundary_names,
        operator_core_dirichlet_names,
        fvm_matrix_names
        ):
    constructor_args = [
        'const std::shared_ptr<const nosh::mesh> & _mesh'
        ]
    init_operator_core_edge = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in operator_core_edge_names]
                )
            )
    init_operator_core_vertex = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in operator_core_vertex_names]
                )
            )
    init_operator_core_boundary = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in operator_core_boundary_names]
                )
            )
    init_operator_core_dirichlet = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>()' % n
                 for n in operator_core_dirichlet_names]
                )
            )
    init_fvm_matrix = '{%s}' % (
            ', '.join(
                ['std::make_shared<%s>(_mesh)' % n
                 for n in fvm_matrix_names]
                )
            )

    members_init = [
      '%s(\n_mesh,\n %s,\n %s,\n %s,\n %s,\n %s\n)' %
      (base_class_name,
       init_operator_core_edge,
       init_operator_core_vertex,
       init_operator_core_boundary,
       init_operator_core_dirichlet,
       init_fvm_matrix
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
