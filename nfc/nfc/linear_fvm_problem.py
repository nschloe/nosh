# -*- coding: utf-8 -*-
#
from .fvm_matrix import gather_dependencies


class LinearFvmProblemCode(object):
    def __init__(self, obj):
        self.class_name = obj.__name__.lower()
        self.dependencies = gather_dependencies(obj)
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_object(self, dep_class_objects):
        code = ''  # TODO
        # code = get_code_linear_problem(
        #     'linear_fvm_problem.tpl',
        #     class_name,
        #     'nosh::linear_problem',
        #     matrix_core_names['edge'],
        #     matrix_core_names['vertex'],
        #     matrix_core_names['boundary'],
        #     matrix_core_names['dirichlet'],
        #     )

        return {
            'code': code
            }
