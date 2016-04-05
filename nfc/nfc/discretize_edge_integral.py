# -*- coding: utf-8 -*-
#
import logging
import sympy
from sympy.matrices.expressions.matexpr import \
        MatrixElement, MatrixExpr, MatrixSymbol


logging.basicConfig(level=logging.DEBUG)


class DiscretizeEdgeIntegral(object):
    def __init__(self):  # , operators, vector_args, scalar_args):
        self.arg_translate = {}
        self._discretization = 0
        self._required_operators = []
        self.x0 = sympy.Symbol('x0')
        self.x1 = sympy.Symbol('x1')
        self.u0 = sympy.Symbol('u0')
        self.u1 = sympy.Symbol('u1')
        self.edge_length = sympy.Symbol('edge_length')
        self.covolume = sympy.Symbol('covolume')
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
                return self.visit_ChainOp(node, sympy.Mul)
            elif node.is_Number:
                return node
            elif node.is_Symbol:
                return node
            elif node.is_Function:
                return self.visit_Call(node)
            elif isinstance(node, MatrixExpr):
                return node

        raise RuntimeError('Unknown node type \"', type(node), '\".')
        return

    def generate(self, node, u, x):
        '''Entrance point to this class.
        '''
        self._discretization = 0
        self._required_operators = []
        out = self.visit(node)
        # Now multiply the value we got out by the covolume
        out *= self.covolume
        # Replace u(x0) by u0, u(x1) by u1.
        out = out.subs(u(self.x0), self.u0)
        out = out.subs(u(self.x1), self.u1)
        # Replace u(x) by 0.5*(u0 + u1) (the edge midpoint)
        out = out.subs(u(x), 0.5 * (self.u0 + self.u1))
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
        print(node)
        print(node.is_Function)
        try:
            id = node.func.__name__
        except AttributeError:
            id = repr(node)
        logging.debug('> Call %s' % id)
        # Handle special functions
        if id == 'dot':
            assert(len(node.args) == 2)
            print(node.args[1])
            print(type(node.args[1]))
            assert(isinstance(node.args[0], MatrixExpr))
            assert(isinstance(node.args[1], MatrixExpr))
            arg0 = self.visit(node.args[0])
            arg1 = self.visit(node.args[1])
            out = node.func(arg0, arg1)
        elif id == 'n_dot_grad':
            assert(len(node.args) == 2)
            f = node.args[0]
            assert(f.is_Function)
            assert(isinstance(node.args[1], MatrixSymbol))
            out = (f(self.x1) - f(self.x0)) / self.edge_length
        else:
            # Default function handling: Assume one argument, e.g., A(x).
            assert(len(node.args) == 1)
            arg = self.visit(node.args[0])
            out = node.func(arg)
        logging.debug('  Call >')
        return out

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
        ret = operator(args[0], args[1])
        for k in range(2, len(args)):
            ret = operator(ret, args[k])
        return ret

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

