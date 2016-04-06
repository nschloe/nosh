# -*- coding: utf-8 -*-
#
import nfl
import os
from string import Template
import sympy
from .helpers import extract_c_expression, templates_dir
from .discretize_edge_integral import DiscretizeEdgeIntegral
from .code_generator_eigen import CodeGeneratorEigen


def is_affine_linear(expr, vars):
    for var in vars:
        if not sympy.Eq(sympy.diff(expr, var, var), 0):
            return False
    return True


# We still need this for pure matrices
# def is_linear(expr, vars):
#     if not is_affine_linear(expr, vars):
#         return False
#     # Check that expr is not affine.
#     if isinstance(expr, int) or isinstance(expr, float):
#         return expr == 0
#     else:
#         return expr.subs([(var, 0) for var in vars]) == 0


class CodeLinearFvmProblem(object):
    def __init__(self, namespace, obj, name):
        self.namespace = namespace
        self.obj = obj
        self.name = name

        if getattr(obj, 'dirichlet_boundary_conditions', None):
            self.dbcs = obj.dirichlet_boundary_conditions
            for dbc in self.dbcs:
                assert(isinstance(dbc, nfl.DirichletBC))
        else:
            self.dbcs = []

        self.dependencies = self.dbcs
        return

    def get_dependencies(self):
        return self.dependencies

    def get_code(self):
        u = sympy.Function('u')
        res = self.obj.eval(u)
        assert(isinstance(res, nfl.Core))
        v, v_affine = self.get_expr_vertex(u, res.vertex)
        edge_expressions, edge_affine = self.get_expr_edge(u, res.edge)

        # check out used an unused symbols for the edge code
        edge_arguments = set([
            sympy.Symbol('x0'),
            sympy.Symbol('x1'),
            sympy.Symbol('edge_length'),
            sympy.Symbol('edge_covolume')
            ])
        edge_unused_arguments = edge_arguments.copy()
        edge_body = []
        edge_used_symbols = set([])
        for e in [edge_expressions[0][0], edge_expressions[0][1],
                  edge_expressions[1][0], edge_expressions[1][1]]:
            edge_used_symbols = edge_used_symbols.union(e.free_symbols)
        edge_unused_arguments -= edge_used_symbols

        edge_undefined_symbols = edge_used_symbols - edge_arguments

        if nfl.n in edge_undefined_symbols or \
           nfl.neg_n in edge_undefined_symbols:
            edge_body.append('const auto n = (x1 - x0) / edge_length;')
            edge_unused_arguments -= set([
                sympy.Symbol('x0'),
                sympy.Symbol('x1')
                ])
            edge_undefined_symbols.remove(nfl.n)
            edge_undefined_symbols.remove(nfl.neg_n)

        assert(len(edge_undefined_symbols) == 0)

        for name in edge_unused_arguments:
            edge_body.append('(void) %s;' % name)

        edge_used_expressions = set()
        for e in [edge_expressions[0][0], edge_expressions[0][1],
                  edge_expressions[1][0], edge_expressions[1][1]]:
            edge_used_expressions = edge_used_expressions.union(
                    [type(atom) for atom in e.atoms(nfl.Expression)]
                    )

        # handle vertex arguments
        vertex_arguments = set([
            sympy.MatrixSymbol('x', 3, 1),
            sympy.Symbol('control_volume')
            ])
        vertex_used_symbols = v.free_symbols
        vertex_unused_arguments = \
            vertex_arguments - vertex_used_symbols.union(v_affine.free_symbols)
        vertex_undefined_symbols = \
            vertex_used_symbols.union(v_affine.free_symbols) - vertex_arguments
        assert(len(vertex_undefined_symbols) == 0)

        vertex_used_expressions = [
                type(atom) for atom in v.atoms(nfl.Expression)
                ]

        members_init = []
        members_declare = []
        used_expressions = edge_used_expressions.union(vertex_used_expressions)
        for expr in used_expressions:
            members_init.append('%s(%s::%s())' % (expr, self.namespace, expr))
            members_declare.append(
                    'const %s::%s %s;' % (self.namespace, expr, expr)
                    )

        if members_init:
            members_init_code = ':\n' + ',\n'.join(members_init)
        else:
            members_init_code = ''

        # template substitution
        with open(os.path.join(templates_dir, 'matrix_core.tpl'), 'r') as f:
            src = Template(f.read())
            matrix_core_code = src.substitute({
                'name': self.name.lower() + '_core',
                'edge00': self.expr_to_code(edge_expressions[0][0]),
                'edge01': self.expr_to_code(edge_expressions[0][1]),
                'edge10': self.expr_to_code(edge_expressions[1][0]),
                'edge11': self.expr_to_code(edge_expressions[1][1]),
                'edge_affine0': self.expr_to_code(-edge_affine[0]),
                'edge_affine1': self.expr_to_code(-edge_affine[1]),
                'edge_body': '\n'.join(edge_body),
                'vertex_contrib': extract_c_expression(v),
                'vertex_affine': extract_c_expression(-v_affine),
                'vertex_body': '\n'.join(
                    ('(void) %s;' % name) for name in vertex_unused_arguments
                    ),
                'members_init': members_init_code,
                'members_declare': '\n'.join(members_declare)
                })

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
            matrix_code = src.substitute({
                'name': self.name.lower(),
                'constructor_args': ',\n'.join(constructor_args),
                'members_init': ',\n'.join(members_init),
                'members_declare': '\n'.join(members_declare)
                })

        code = matrix_core_code + '\n' + matrix_code
        return code

    def expr_to_code(self, expr):
        gen = CodeGeneratorEigen()
        code, _ = gen.generate(expr)
        return code

    def get_expr_vertex(self, u, function):
        # Numerically integrate function over a control volume.
        # TODO find out if the code can be evaluated at the volume "midpoint",
        # then evaluate it there, and multiply by the volume.
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
        # Make sure that f is affine linear in u0; we're building a matrix
        # here.
        assert(is_affine_linear(fu0, [u0]))
        control_volume = sympy.Symbol('control_volume')
        # Get coefficient of u0
        coeff = sympy.diff(fu0, u0)
        coeff = control_volume * coeff
        # Get affine part
        if isinstance(fu0, float):
            affine = fu0
        else:
            affine = fu0.subs(u0, 0)
        return (coeff, affine)

    def get_expr_edge(self, u, function):
        x = sympy.MatrixSymbol('x', 3, 1)

        generator = DiscretizeEdgeIntegral()
        expr, _ = generator.generate(function(x), u, x)

        # Make sure the expression is linear in u0, u1.
        if not is_affine_linear(expr, [generator.u0, generator.u1]):
            raise RuntimeError((
                'The given function\n'
                '    f(x) = %s\n'
                'does not seem to be affine linear in u.')
                % function(x)
                )

        # Get the coefficients of u0, u1.
        coeff00 = sympy.diff(expr, generator.u0)
        coeff01 = sympy.diff(expr, generator.u1)

        # Now construct the coefficients for the other way around.
        coeff10 = coeff01.subs([
            (generator.u0, generator.u1),
            (generator.u1, generator.u0),
            (nfl.n, nfl.neg_n)
            ])
        coeff11 = coeff00.subs([
            (generator.u0, generator.u1),
            (generator.u1, generator.u0),
            (nfl.n, nfl.neg_n)
            ])

        affine = expr.subs([(generator.u0, 0), (generator.u1, 0)])

        return (
            [[coeff00, coeff01], [coeff10, coeff11]],
            [affine, affine]
            )
