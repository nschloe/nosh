# -*- coding: utf-8 -*-
#
class CodeMatrixCore(object):
    def __init__(self, namespace, core, name):
        # handle the edge contributions
        if callable(getattr(core, 'edge_contrib', None)):
            edge_result, edge_unused_symbols, edge_used_expressions = \
                    self.get_expr_edge_contrib(core.edge_contrib)

        if callable(getattr(core, 'vertex_contrib', None)):
            vertex_result, vertex_unused_symbols, vertex_used_expressions = \
                    self.get_expr_vertex_contrib(core.vertex_contrib)

        # Go through all used expressions and check if they are defined. If
        # not, they need to be given as arguments to the constructor.
        members_init = []
        members_declare = []
        used_expressions = edge_used_expressions.union(vertex_used_expressions)
        for expr in used_expressions:
            members_init.append('%s(%s::%s())' % (expr, namespace, expr))
            members_declare.append(
                    'const %s::%s %s;' % (namespace, expr, expr)
                    )

        if members_init:
            members_init_code = ':\n' + ',\n'.join(members_init)
        else:
            members_init_code = ''

        self.dependencies = used_expressions

        # template substitution
        with open(os.path.join(templates_dir, 'matrix_core.tpl'), 'r') as f:
            src = Template(f.read())
            self.code = src.substitute({
                'name': name.lower(),  # class names are lowercase
                'edge00': extract_c_expression(edge_result[0][0]),
                'edge01': extract_c_expression(edge_result[0][1]),
                'edge10': extract_c_expression(edge_result[1][0]),
                'edge11': extract_c_expression(edge_result[1][1]),
                'edge_body': '\n'.join(
                    ('(void) %s;' % name) for name in edge_unused_symbols
                    ),
                'vertex_contrib': extract_c_expression(vertex_result),
                'vertex_body': '\n'.join(
                    ('(void) %s;' % name) for name in vertex_unused_symbols
                    ),
                'members_init': members_init_code,
                'members_declare': '\n'.join(members_declare)
                })

    def get_dependencies(self):
        return self.dependencies

    def get_code(self):
        return self.code

    def get_expr_edge_contrib(self, method):
        x0 = sympy.MatrixSymbol('x0', 3, 1)
        x1 = sympy.MatrixSymbol('x1', 3, 1)
        edge_length = sympy.Symbol('edge_length')
        edge_covolume = sympy.Symbol('edge_covolume')
        all_symbols = set([x0, x1, edge_length, edge_covolume])

        specs = inspect.getargspec(method)
        assert(len(specs.args) == len(all_symbols) + 1)

        result = method(x0, x1, edge_length, edge_covolume)

        assert(len(result) == 2)
        assert(len(result[0]) == 2)
        assert(len(result[1]) == 2)

        # Check if any of the arguments is not used in the function.
        # (We'll declare them (void) to supress compiler warnings.)
        used_symbols = set()
        used_expressions = set()
        for r in [result[0][0], result[0][1], result[1][0], result[1][1]]:
            used_expressions = used_expressions.union(
                    [type(atom) for atom in r.atoms(nfl.Expression)]
                    )
            used_symbols = used_symbols.union(r.free_symbols)

        unused_args = all_symbols - used_symbols
        return result, unused_args, used_expressions

    def get_expr_vertex_contrib(self, method):
        # handle the vertex contributions
        x = sympy.MatrixSymbol('x')
        vol = sympy.Symbol('control_volume')
        all_symbols = set([x, vol])

        specs = inspect.getargspec(method)
        assert(len(specs.args) == len(all_symbols) + 1)

        result = method(x, vol)

        try:
            unused_args = all_symbols - result.free_symbols
        except AttributeError:
            unused_args = all_symbols

        try:
            used_expressions = set(
                    [type(atom) for atom in result.atoms(nfl.Expression)]
                    )
        except AttributeError:
            used_expressions = set([])

        return result, unused_args, used_expressions


