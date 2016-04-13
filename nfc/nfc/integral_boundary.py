# -*- coding: utf-8 -*-
#
import os
from string import Template
import sympy
import nfl

from .code_generator_eigen import CodeGeneratorEigen
from .helpers import \
        extract_c_expression, \
        extract_linear_components, \
        get_uuid, \
        is_affine_linear, \
        list_unique, \
        members_init_declare, \
        templates_dir


class IntegralBoundary(object):
    def __init__(self, namespace, u, integrand, subdomains, is_matrix):
        self.namespace = namespace

        self.is_matrix = is_matrix

        self.class_name = 'matrix_core_boundary_' + get_uuid()

        self.expr, self.u0 = \
            _discretize_integral(u, integrand)

        self.dependencies = set().union(
            [type(atom) for atom in self.expr.atoms(nfl.Expression)],
            subdomains
            )
        return

    def get_dependencies(self):
        return self.dependencies

    def get_class_object(self, dep_class_objects):
        arguments = set([sympy.Symbol('vertex')])
        used_vars = self.expr.free_symbols
        if self.u0 in used_vars:
            used_vars.remove(self.u0)
        extra_body, extra_init, extra_declare = \
                _get_extra(arguments, used_vars)

        # now take care of the template substitution
        members_init, members_declare = \
            members_init_declare(
                    self.namespace,
                    'matrix_core_boundary',
                    dep_class_objects
                    )

        members_init.extend(extra_init)
        members_declare.extend(extra_declare)

        if members_init:
            members_init_code = ':\n' + ',\n'.join(members_init)
        else:
            members_init_code = ''

        if self.is_matrix:
            coeff, affine = extract_linear_components(self.expr, self.u0)
            type = 'matrix_core_boundary'
            filename = os.path.join(templates_dir, 'matrix_core_boundary.tpl')
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'coeff': extract_c_expression(coeff),
                    'affine': extract_c_expression(-affine),
                    'body': '\n'.join(extra_body),
                    'members_init': members_init_code,
                    'members_declare': '\n'.join(members_declare)
                    })
        else:
            type = 'matrix_core_operator'
            filename = os.path.join(
                    templates_dir, 'operator_core_boundary.tpl'
                    )
            with open(filename, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': self.class_name,
                    'coeff': extract_c_expression(boundary_coeff),
                    'affine': extract_c_expression(-boundary_affine),
                    'body': '\n'.join(
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


# def get_matrix_core_boundary_code(namespace, class_name, core):
#     '''Get code generator from raw core object.
#     '''
#     # handle the boundary contributions
#     x = sympy.MatrixSymbol('x')
#     vol = sympy.Symbol('control_volume')
#     all_symbols = set([x, vol])
#
#     specs = inspect.getargspec(method)
#     assert(len(specs.args) == len(all_symbols) + 1)
#
#     boundary_coeff, boundary_affine = method(x, vol)
#
#     return _get_code_matrix_core_boundary(
#             namespace, class_name,
#             boundary_coeff, boundary_affine
#             )


def _discretize_integral(u, function):
    # Numerically integrate function over a control volume.
    x = sympy.MatrixSymbol('x', 3, 1)
    # Evaluate the function for u at x.
    fx = function(x)
    # Replace all occurences of u(x) by u0 (the value at the control volume
    # center) and multiply by the control volume)
    u0 = sympy.Symbol('u0')
    try:
        fu0 = fx.subs(u(x), u0)
    except AttributeError:
        # 'int' object has no attribute 'subs'
        fu0 = fx
    surface_area = sympy.Symbol('surface_area')
    return surface_area * fu0, u0


def _get_extra(arguments, used_variables):
    vertex = sympy.Symbol('vertex')
    unused_arguments = arguments - used_variables
    undefined_symbols = used_variables - arguments

    init = []
    body = []
    declare = []

    surface_area = sympy.Symbol('surface_area')
    if surface_area in undefined_symbols:
        init.append('mesh_(mesh)')
        declare.append('const std::shared_ptr<const nosh::mesh> mesh_;')
        init.append('surfs_(mesh->boundary_surface_areas())')
        declare.append('const std::vector<double> surfs_;')
        body.append('const auto k = this->mesh_->local_index(vertex);')
        body.append('const auto surface_area = this->surfs_[k];')
        undefined_symbols.remove(surface_area)
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
