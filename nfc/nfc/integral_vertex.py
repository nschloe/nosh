# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .code_generator_eigen import CodeGeneratorEigen
from .expression import *
from .subdomain import *
from .helpers import \
        extract_c_expression, \
        extract_linear_components, \
        get_uuid, \
        is_affine_linear, \
        list_unique, \
        members_init_declare, \
        templates_dir


class IntegralVertex(object):
    def __init__(self, namespace, integrand, subdomains, matrix_var=None):
        self.namespace = namespace
        self.class_name = 'vertex_core_' + get_uuid()

        self.matrix_var = matrix_var

        x = sympy.MatrixSymbol('x', 3, 1)
        fx = integrand(x)
        self.expr, self.vector_vars = _discretize_expression(fx)

        self.dependencies = set().union(
            [ExpressionCode(type(atom))
                for atom in self.expr.atoms(nfl.Expression)],
            [SubdomainCode(sd) for sd in subdomains]
            )
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_object(self, dependency_class_objects):
        if self.matrix_var:
            parent_class = 'matrix_core_vertex'
        else:
            parent_class = 'operator_core_vertex'

        arguments = set([sympy.Symbol('vertex')])
        used_vars = self.expr.free_symbols
        for vector_var in self.vector_vars:
            if vector_var in used_vars:
                used_vars.remove(vector_var)
        extra_body, extra_init, extra_declare = _get_extra(
                arguments, used_vars
                )

        # now take care of the template substitution
        members_init, members_declare = \
            members_init_declare(
                    self.namespace,
                    parent_class,
                    dependency_class_objects
                    )

        members_init.extend(extra_init)
        members_declare.extend(extra_declare)

        if members_init:
            members_init_code = ':\n' + ',\n'.join(members_init)
        else:
            members_init_code = ''

        if self.matrix_var:
            coeff, affine = extract_linear_components(
                    self.expr,
                    sympy.Symbol('%s[k]' % self.matrix_var)
                    )
            type = 'matrix_core_vertex'
            filename = os.path.join(templates_dir, 'matrix_core_vertex.tpl')
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'vertex_contrib': extract_c_expression(coeff),
                    'vertex_affine': extract_c_expression(-affine),
                    'vertex_body': '\n'.join(extra_body),
                    'members_init': members_init_code,
                    'members_declare': '\n'.join(members_declare)
                    })
        else:
            type = 'operator_core_vertex'
            filename = os.path.join(templates_dir, 'operator_core_vertex.tpl')
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'return_value': extract_c_expression(self.expr),
                    'eval_body': '\n'.join(
                        ('(void) %s;' % name) for name in unused_args
                        ),
                    'members_init': members_init_code,
                    'members_declare': '\n'.join(members_declare)
                    })

        return {
            'type': type,
            'code': code,
            'class_name': self.class_name,
            'constructor_args': []
            }


# def get_matrix_core_vertex_code(namespace, class_name, core):
#     '''Get code generator from raw core object.
#     '''
#     # handle the vertex contributions
#     x = sympy.MatrixSymbol('x')
#     vol = sympy.Symbol('control_volume')
#     all_symbols = set([x, vol])
#
#     specs = inspect.getargspec(method)
#     assert(len(specs.args) == len(all_symbols) + 1)
#
#     vertex_coeff, vertex_affine = method(x, vol)
#
#     return _get_code_matrix_core_vertex(
#             namespace, class_name,
#             vertex_coeff, vertex_affine
#             )


def _discretize_expression(expr):
    # Find all Nosh function variables
    fks = []
    if isinstance(expr, float) or isinstance(expr, int):
        pass
    else:
        function_vars = []
        for f in expr.atoms(sympy.Function):
            if hasattr(f, 'nosh'):
                function_vars.append(f)

        for function_var in function_vars:
            # Replace all occurences of u(x) by u[k] (the value at the control
            # volume center) and multiply by the control volume)
            fk = sympy.Symbol('%s[k]' % function_var.func)
            try:
                expr = expr.subs(function_var, fk)
            except AttributeError:
                # 'int' object has no attribute 'subs'
                pass
            fks.append(fk)

    return sympy.Symbol('control_volume') * expr, fks


def _get_extra(arguments, used_variables):
    vertex = sympy.Symbol('vertex')
    unused_arguments = arguments - used_variables
    undefined_symbols = used_variables - arguments

    init = []
    body = []
    declare = []

    control_volume = sympy.Symbol('control_volume')
    if control_volume in undefined_symbols:
        init.append('mesh_(mesh)')
        declare.append('const std::shared_ptr<const nosh::mesh> mesh_;')
        init.append('c_data_(mesh->control_volumes()->getData())')
        declare.append('const Teuchos::ArrayRCP<const double> c_data_;')
        body.append('const auto k = this->mesh_->local_index(vertex);')
        body.append('const auto control_volume = this->c_data_[k];')
        undefined_symbols.remove(control_volume)
        if vertex in unused_arguments:
            unused_arguments.remove(vertex)

    x = sympy.MatrixSymbol('x', 3, 1)
    if x in undefined_symbols:
        init.append('mesh_(mesh)')
        declare.append('const std::shared_ptr<const nosh::mesh> mesh_;')
        init.append('c_data_(mesh->control_volumes()->getData())')
        declare.append('const Teuchos::ArrayRCP<const double> c_data_;')
        body.append('const auto k = this->mesh_->local_index(vertex);')
        body.append('const auto x = this->mesh_->get_coords(vertex);')
        undefined_symbols.remove(x)
        if vertex in unused_arguments:
            unused_arguments.remove(vertex)

    if len(undefined_symbols) > 0:
        raise RuntimeError(
                'The following symbols are undefined: %s' % undefined_symbols
                )

    # remove double lines
    body = list_unique(body)
    init = list_unique(init)
    declare = list_unique(declare)

    for name in unused_arguments:
        body.insert(0, '(void) %s;' % name)

    return body, init, declare
