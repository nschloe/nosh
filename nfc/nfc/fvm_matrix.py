# -*- coding: utf-8 -*-
#
import os
import sympy


def get_code_fvm_matrix(name, obj):
    '''Get Nosh C++ code for an FVM matrix.
    '''
    if getattr(obj, 'matrix_cores', None):
        matrix_cores = obj.matrix_cores
        for matrix_core in matrix_cores:
            assert(isinstance(matrix_core, nfl.MatrixCore))
    else:
        matrix_cores = []

    if getattr(obj, 'boundary_conditions', None):
        bcs = obj.boundary_conditions
        for bc in bcs:
            assert(isinstance(bc, nfl.DirichletBC))
    else:
        boundary_conditions = []

    dependencies = matrix_cores + bcs

    return _get_code(name, matrix_cores, bcs), dependencies


def _get_code(name, matrix_cores, bcs):
    # handle the matrix cores
    joined = ', '.join(
            'std::make_shared<%s>()' %
            type(mc).__name__.lower() for mc in matrix_cores
            )
    code_matrix_cores = '{%s}' % joined
    # matrix cores handling done

    # handle the boundary conditions
    joined = ', '.join(
            'std::make_shared<%s>()' %
            type(bc).__name__.lower() for bc in bcs
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
            'name': name.lower(),  # class names are lowercase
            'constructor_args': ',\n'.join(constructor_args),
            'members_init': ',\n'.join(members_init),
            'members_declare': '\n'.join(members_declare)
            })

    return code
