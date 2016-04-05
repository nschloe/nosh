# -*- coding: utf-8 -*-
#
import logging
import sympy
from sympy.matrices.expressions.matexpr import MatrixElement


logging.basicConfig(level=logging.DEBUG)


class DiscretizeEdgeIntegral(object):
    def __init__(self):  # , operators, vector_args, scalar_args):
        self.arg_translate = {}
        self._discretization = 0
        self._required_operators = []
        return

    def visit(self, node):
        if isinstance(node, int):
            return node
        elif isinstance(node, float):
            return node
        elif isinstance(node, sympy.Basic):
            if node.is_Add:
                return self.visit_ChainOp(node, sympy.Add)
            elif node.is_Mul:
                return self.visit_ChainOp(node, sympy.Multiply)
            elif node.is_Number:
                return node
            elif node.is_Symbol:
                return node
            elif node.is_Function:
                a = self.visit_Call(node)
                return a
            elif isinstance(node, MatrixElement):
                return node

        raise RuntimeError('Unknown node type \"', type(node), '\".')
        return

    def generate(self, node):
        '''Entrance point to this class.
        '''
        self._discretization = 0
        self._required_operators = []
        out = self.visit(node)
        print('out', out)
        exit()
        return out, self._required_operators

    def generic_visit(self, node):
        raise RuntimeError(
            'Should never be called. __name__:', type(node).__name__
            )
        self.visit(node)
        return

    def visit_Load(self, node):
        logging.debug('> Load >')
        pass

    def visit_Call(self, node):
        '''Handles calls for operators A(u) and pointwise functions sin(u).
        '''
        id = node.func.__name__
        logging.debug('> Call %s' % id)
        # Check if this is the top (or root) of the recursion. If it is, the
        # output variable will be `y`.
        assert(len(node.args) == 1)  # one argument, e.g., A(x)
        ret = self.visit(node.args[0])
        if isinstance(node, nfl.FvmMatrix):
            # The argument to A(.) must be of vector type.
            if isinstance(ret, Vector):
                arg_name = ret
            elif isinstance(ret, Pointwise):
                arg_name = self._get_outvector()
                self._to_vector(ret, arg_name)
            else:
                raise ValueError('Illegal input type')
            # Get the output vector
            # Put it all together
            self._code += '\n%s->apply(%s, %s);\n' \
                % (id + '_', arg_name, out_vector)
            self._required_operators.append({
                'var_name': id + '_',
                'class_name': id,
                })
            logging.debug('  Call >')
            return Vector(out_vector)
        else:
            # Assume that the operator is a C++ intrinsic or otherwise defined.
            # The argument must be of pointwise type.
            if isinstance(ret, Vector):
                a = self._to_pointwise(ret)
            elif isinstance(ret, Pointwise):
                a = ret
            else:
                raise ValueError('Illegal input type')
            logging.debug('  Call >')
            return Pointwise('%s(%s)' % (id, a))

    def visit_UnaryOp(self, node):
        '''Handles unary operations (e.g., +, -,...).
        '''
        code_op = self.visit(node.op)
        logging.debug('> UnaryOp %s' % code_op)
        ret = self.visit(node.operand)
        if isinstance(ret, Vector):
            code = self._to_pointwise(ret)
        elif isinstance(ret, Pointwise):
            code = ret
        else:
            raise ValueError('Illegal input type')
        # plug it together
        pointwise_code = '%s%s' % (code_op, code)
        logging.debug('  UnaryOp >')
        return Pointwise(pointwise_code)

    def visit_ChainOp(self, node, operator):
        '''Handles binary operations (e.g., +, -, *,...).
        '''
        logging.debug('> BinOp %s' % operator)
        # collect the pointwise code for left and right
        args = []
        for n in node.args:
            ret = self.visit(n)
            args.append(ret)
        # plug it together
        return operator(args[0], args[1])

    def visit_Name(self, node):
        logging.debug('> Name %s >' % id)
        return node.name

    def visit_Add(self, node):
        return '+'

    def visit_Sub(self, node):
        return '-'

    def visit_Mult(self, node):
        return '*'

    def visit_Div(self, node):
        return '/'

    def visit_UAdd(self, node):
        return '+'

    def visit_USub(self, node):
        return '-'

    def visit_Num(self, node):
        return Pointwise(str(node.n))

