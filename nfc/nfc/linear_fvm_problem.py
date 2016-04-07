# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
from .matrix_core import get_core_code_from_expression
from .helpers import templates_dir


class CodeLinearFvmProblem(object):
    def __init__(self, namespace, name, obj):
        self.name = name

        if getattr(obj, 'dirichlet_boundary_conditions', None):
            self.dbcs = set(obj.dirichlet_boundary_conditions)
            for dbc in self.dbcs:
                assert(isinstance(dbc, nfl.DirichletBC))
        else:
            self.dbcs = set()

        self.dependencies = self.dbcs

        core_code_gen = get_core_code_from_expression(
                namespace, name, obj.eval
                )
        self.dependencies.update(core_code_gen.get_dependencies())
        self.code = core_code_gen.get_code() + \
            '\n' + \
            self.get_code_linear_problem()
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
        init_matrix_cores = '{std::make_shared<%s>()}' % (
                self.name.lower() + '_core'
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
