# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
from .helpers import \
        compare_variables, \
        extract_c_expression, \
        sanitize_identifier, \
        templates_dir


class ExpressionCode(object):
    def __init__(self, cls):
        self.class_name = sanitize_identifier(cls.__name__)
        self.cls = cls
        return

    def get_dependencies(self):
        return []

    def get_class_object(self, dependency_class_objects):
        x = sympy.MatrixSymbol('x', 3, 1)
        result = self.cls.eval(x)
        if result is None:
            # The expression must still be defined
            code = ''
        else:
            unused_args, _ = compare_variables(set([x]), [result])

            # template substitution
            filename = os.path.join(templates_dir, 'expression.tpl')
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'eval': extract_c_expression(result),
                    'unused_args': '\n'.join(
                        ('(void) %s;' % name) for name in unused_args
                        ),
                    })

        return {
            'code': code,
            'type': 'expression',
            'class_name': self.class_name
            }
