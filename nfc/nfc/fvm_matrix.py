# -*- coding: utf-8 -*-
#
import os
import sympy


class CodeFvmMatrix(object):
    '''Get Nosh C++ code for an FVM matrix.
    '''
    def __init__(self, eo, name):
        self.eo = eo
        self.name = name

        if getattr(eo, 'matrix_cores', None):
            self.matrix_cores = eo.matrix_cores
            for matrix_core in self.matrix_cores:
                assert(isinstance(matrix_core, nfl.MatrixCore))
        else:
            self.matrix_cores = []

        if getattr(eo, 'boundary_conditions', None):
            self.bcs = eo.boundary_conditions
            for bc in self.bcs:
                assert(isinstance(bc, nfl.DirichletBC))
        else:
            self.boundary_conditions = []

        self.dependencies = self.matrix_cores + self.bcs

    def get_dependencies(self):
        return self.dependencies

    def get_code(self):
        # handle the matrix cores
        joined = ', '.join(
                'std::make_shared<%s>()' %
                type(mc).__name__.lower() for mc in self.matrix_cores
                )
        code_matrix_cores = '{%s}' % joined
        # matrix cores handling done

        # handle the boundary conditions
        joined = ', '.join(
                'std::make_shared<%s>()' %
                type(bc).__name__.lower() for bc in self.bcs
                )
        code_bcs = '{%s}' % joined
        # boundary conditions handling done

        constructor_args = [
            'const std::shared_ptr<const nosh::mesh> & _mesh'
            ]
        members_init = [
          'nosh::fvm_matrix(\n_mesh,\n %s,\n %s\n)' %
          (code_matrix_cores, code_bcs)
          ]
        members_declare = []

        # template substitution
        with open(os.path.join(templates_dir, 'fvm_matrix.tpl'), 'r') as f:
            src = Template(f.read())
            code = src.substitute({
                'name': self.name.lower(),  # class names are lowercase
                'constructor_args': ',\n'.join(constructor_args),
                'members_init': ',\n'.join(members_init),
                'members_declare': '\n'.join(members_declare)
                })

        return code


