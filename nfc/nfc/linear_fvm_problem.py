# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
import sympy

from .edge_core import get_edge_core_code_from_expression
from .vertex_core import get_vertex_core_code_from_expression
from .boundary_core import get_boundary_core_code_from_expression
from .helpers import templates_dir


class CodeLinearFvmProblem(object):
    def __init__(self, namespace, name, obj):
        self.name = name
        self.code = ''

        if getattr(obj, 'dirichlet_boundary_conditions', None):
            self.dbcs = set(obj.dirichlet_boundary_conditions)
            for dbc in self.dbcs:
                assert(isinstance(dbc, nfl.DirichletBC))
        else:
            self.dbcs = set()

        self.dependencies = self.dbcs

        u = sympy.Function('u')
        res = obj.eval(u)
        assert(isinstance(res, nfl.Core))
        for core in res.cores:
            if core[1] == 'dS':
                core_code_gen = get_edge_core_code_from_expression(
                        namespace, name, u, core[0]
                        )
            elif core[1] == 'dV':
                core_code_gen = get_vertex_core_code_from_expression(
                        namespace, name, u, core[0]
                        )
            elif core[1] == 'dGamma':
                core_code_gen = get_boundary_core_code_from_expression(
                        namespace, name, u, core[0]
                        )
            else:
                raise RuntimeError('Illegal core type \'%s\'.' % core[1])

            self.dependencies.update(core_code_gen.get_dependencies())

            self.code += '\n' + core_code_gen.get_code()

        self.code += '\n' + self.get_code_linear_problem()
        return

    def get_dependencies(self):
        return self.dependencies

    def get_code(self):
        return self.code

    def get_code_linear_problem(self):
        # fvm_matrix code
        constructor_args = [
            'const std::shared_ptr<const nosh::mesh> & _mesh'
            ]
        init_edge_cores = '{std::make_shared<%s>()}' % (
                self.name.lower() + '_edge_core'
                )

        # handle the boundary conditions
        joined = ', '.join(
                'std::make_shared<%s>()' %
                type(bc).__name__.lower() for bc in self.dbcs
                )
        init_dbcs = '{%s}' % joined
        # boundary conditions handling done

        members_init = [
          'nosh::linear_problem(\n_mesh,\n %s,\n %s\n)' %
          (init_matrix_cores, init_dbcs)
          ]
        members_declare = []

        templ = os.path.join(templates_dir, 'linear_fvm_problem.tpl')
        with open(templ, 'r') as f:
            src = Template(f.read())
            code = src.substitute({
                'name': self.name.lower(),
                'constructor_args': ',\n'.join(constructor_args),
                'members_init': ',\n'.join(members_init),
                'members_declare': '\n'.join(members_declare)
                })

        return code
