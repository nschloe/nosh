# -*- coding: utf-8 -*-
#
from .dirichlet import *
from .integral_boundary import *
from .integral_edge import *
from .integral_vertex import *
from .helpers import sanitize_identifier
from .fvm_matrix import gather_core_dependencies, get_code_linear_problem


class LinearFvmProblemCode(object):
    def __init__(self, namespace, cls):
        self.class_name = sanitize_identifier(cls.__name__)
        self.namespace = namespace

        u = sympy.Function('u')
        res = cls.apply(u)
        self.dependencies = \
            gather_core_dependencies(
                    namespace, res, cls.dirichlet, is_matrix=True
                    )
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_object(self, dep_class_objects):

        # Go through the dependencies collect the cores.
        dirichlet_core_names = []
        vertex_core_names = []
        edge_core_names = []
        boundary_core_names = []
        for dep in self.dependencies:
            if isinstance(dep, Dirichlet):
                dirichlet_core_names.append(dep.class_name)
            elif isinstance(dep, IntegralVertex):
                vertex_core_names.append(dep.class_name)
            elif isinstance(dep, IntegralEdge):
                edge_core_names.append(dep.class_name)
            elif isinstance(dep, IntegralBoundary):
                boundary_core_names.append(dep.class_name)

        code = get_code_linear_problem(
            'linear_fvm_problem.tpl',
            self.class_name,
            'nosh::linear_problem',
            edge_core_names,
            vertex_core_names,
            boundary_core_names,
            dirichlet_core_names
            )

        return {
            'code': code
            }
