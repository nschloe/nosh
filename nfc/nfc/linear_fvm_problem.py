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

        # edge
        edge_coeff, edge_affine = self.get_expr_edge(u, res.edge)
        edge_body, edge_used_expressions = \
            self.get_code_elements_edge(edge_coeff, edge_affine)

        # vertex
        vertex_coeff, vertex_affine = self.get_expr_vertex(u, res.vertex)
        vertex_arguments = set([
            sympy.MatrixSymbol('x', 3, 1),
            sympy.Symbol('control_volume')
            ])
        vertex_unused_arguments, vertex_used_expressions = \
            self.scan_code(vertex_arguments, [vertex_coeff, vertex_affine])

        # domain boundary
        db_coeff, db_affine = self.get_expr_db(res.domain_boundary)
        db_arguments = set([
            sympy.MatrixSymbol('x', 3, 1),
            sympy.Symbol('surface_area')
            ])
        db_unused_arguments, db_used_expressions = \
            self.scan_code(db_arguments, [db_coeff, db_affine])

        # now take care of the template substitution
        members_init = []
        members_declare = []

        # init and declare all expressions
        used_expressions = set().union(
                edge_used_expressions,
                vertex_used_expressions,
                db_used_expressions
                )
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
                'edge00': self.expr_to_code(edge_coeff[0][0]),
                'edge01': self.expr_to_code(edge_coeff[0][1]),
                'edge10': self.expr_to_code(edge_coeff[1][0]),
                'edge11': self.expr_to_code(edge_coeff[1][1]),
                'edge_affine0': self.expr_to_code(-edge_affine[0]),
                'edge_affine1': self.expr_to_code(-edge_affine[1]),
                'edge_body': '\n'.join(edge_body),
                'vertex_contrib': extract_c_expression(vertex_coeff),
                'vertex_affine': extract_c_expression(-vertex_affine),
                'vertex_body': '\n'.join(
                    ('(void) %s;' % name) for name in vertex_unused_arguments
                    ),
                'db_coeff': extract_c_expression(db_coeff),
                'db_affine': extract_c_expression(-db_affine),
                'db_body': '\n'.join(
                    ('(void) %s;' % name) for name in db_unused_arguments
                    ),
                'members_init': members_init_code,
                'members_declare': '\n'.join(members_declare)
                })

        linear_problem_code = self.get_linear_problem_code()

        return matrix_core_code + '\n' + linear_problem_code

    def get_linear_problem_code(self):
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
            affine = control_volume * fu0
        else:
            affine = control_volume * fu0.subs(u0, 0)
        return (coeff, affine)

    def get_expr_db(self, function):
        # Numerically integrate function over a part of the domain boundary,
        # namely the part surrounding a vertex.
        x = sympy.MatrixSymbol('x', 3, 1)
        # Evaluate the function for u at x.
        fx = function(x)
        # Replace all occurences of u(x) by u0 (the value at the center) and
        # multiply by the surface area)
        u0 = sympy.Symbol('u0')
        try:
            fu0 = fx.subs(u(x), u0)
        except AttributeError:
            # 'int' object has no attribute 'subs'
            fu0 = fx
        # Make sure that f is affine linear in u0; we're building a matrix
        # here.
        assert(is_affine_linear(fu0, [u0]))
        control_volume = sympy.Symbol('surface_area')
        # Get coefficient of u0
        coeff = sympy.diff(fu0, u0)
        coeff = control_volume * coeff
        # Get affine part
        if isinstance(fu0, float) or isinstance(fu0, int):
            affine = control_volume * fu0
        else:
            affine = control_volume * fu0.subs(u0, 0)
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

    def get_code_elements_edge(self, edge_coeff, edge_affine):
        # Check out used an unused symbols for the edge code.
        # Arguments from the template.
        arguments = set([
            sympy.Symbol('x0'),
            sympy.Symbol('x1'),
            sympy.Symbol('edge_length'),
            sympy.Symbol('edge_covolume')
            ])
        unused_arguments = arguments.copy()
        body = []
        used_symbols = set([])
        for e in [edge_coeff[0][0], edge_coeff[0][1],
                  edge_coeff[1][0], edge_coeff[1][1],
                  edge_affine[0], edge_affine[1]
                  ]:
            used_symbols = used_symbols.union(e.free_symbols)
        unused_arguments -= used_symbols

        undefined_symbols = used_symbols - arguments

        if nfl.n in undefined_symbols or \
           nfl.neg_n in undefined_symbols:
            body.append('const auto n = (x1 - x0) / edge_length;')
            unused_arguments -= set([
                sympy.Symbol('x0'),
                sympy.Symbol('x1')
                ])
            undefined_symbols.remove(nfl.n)
            undefined_symbols.remove(nfl.neg_n)

        assert(len(undefined_symbols) == 0)

        for name in unused_arguments:
            body.append('(void) %s;' % name)

        used_expressions = set()
        for e in [edge_coeff[0][0], edge_coeff[0][1],
                  edge_coeff[1][0], edge_coeff[1][1],
                  edge_affine[0], edge_affine[1]
                  ]:
            used_expressions = used_expressions.union(
                    [type(atom) for atom in e.atoms(nfl.Expression)]
                    )

        return body, used_expressions

    def scan_code(self, arguments, expressions):
        used_symbols = set([])
        used_expressions = set([])

        for expr in expressions:
            try:
                used_symbols.update(expr.free_symbols)
            except AttributeError:
                pass

            used_expressions.update(set([
                    type(atom) for atom in expr.atoms(nfl.Expression)
                    ]))

        unused_arguments = arguments - used_symbols
        undefined_symbols = used_symbols - arguments
        assert(len(undefined_symbols) == 0)

        return unused_arguments, used_expressions
