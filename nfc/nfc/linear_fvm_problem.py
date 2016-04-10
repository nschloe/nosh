# -*- coding: utf-8 -*-
#
from .fvm_matrix import handle_core_dependencies, get_code_linear_problem


def get_code_linear_fvm_problem(namespace, class_name, obj):

    code, dependencies, matrix_core_names = \
            handle_core_dependencies(namespace, obj)

    # append the code of the linear problem itself
    code += '\n' + get_code_linear_problem(
            'linear_fvm_problem.tpl',
            class_name.lower(),
            'nosh::linear_problem',
            matrix_core_names['edge'],
            matrix_core_names['vertex'],
            matrix_core_names['boundary'],
            matrix_core_names['dirichlet'],
            )

    return code, dependencies
