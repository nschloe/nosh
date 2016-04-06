# -*- coding: utf-8 -*-
#
class CodeOperator(object):
    def __init__(self, operator, name):
        self.operator = operator
        self.name = name
        return

    def get_code(self):
        assert(isinstance(operator, nfl.Operator))

        x = sympy.MatrixSymbol('x', 3, 1)

        assert(callable(operator.eval))

        if (len(inspect.getargspec(operator.eval).args) == 1):
            generator = CodeGen()
            u = sympy.Symbol('u')
            op_code, required_ops = \
                generator.generate(operator.eval(u), {'u': 'x'})

            members = ['const std::shared_ptr<const nosh::mesh> mesh_;']
            members_init = ['mesh_(mesh)']
            for required_operator in required_ops:
                members.append(
                        'const std::shared<const Tpetra::Operator> %s;' %
                        required_operator['var_name'].lower()
                        )
                members_init.append(
                        '%s(std::make_shared<%s>(mesh))' %
                        (required_operator['var_name'].lower(),
                            required_operator['class_name'].lower())
                        )

            # TODO
            # Check if any of the arguments is not used in the function.
            # (We'll declare them (void) to supress compiler warnings.)

            # template substitution
            with open(os.path.join(templates_dir, 'operator.tpl'), 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': name.lower(),  # class names are lowercase
                    'apply': op_code,
                    'members': '\n'.join(members),
                    'members_init': ',\n'.join(members_init)
                    })
        elif (len(inspect.getargspec(operator.eval).args) == 2):
            code = ''
            generator = CodeGen()
            u = sympy.Symbol('u')
            u0 = sympy.Symbol('u0')
            op_code, required_ops = \
                generator.generate(operator.eval(u, u0), {'u': 'x', 'u0': 'x0_'})

            members = [
                    'const std::shared_ptr<const nosh::mesh> mesh_;'
                    'Tpetra::Vector<double,int,int> x0_;'
                    ]
            members_init = ['mesh_(mesh)', 'x0_(x0)']
            for required_operator in required_ops:
                members.append(
                        'const std::shared<const Tpetra::Operator> %s;' %
                        required_operator['var_name'].lower()
                        )
                members_init.append(
                        '%s(std::make_shared<%s>(mesh))' %
                        (required_operator['var_name'].lower(),
                            required_operator['class_name'].lower())
                        )
            # template substitution
            template = os.path.join(templates_dir, 'operator_with_rebuild.tpl')
            with open(template, 'r') as f:
                src = Template(f.read())
                code = src.substitute({
                    'name': name.lower(),  # class names are lowercase
                    'apply': op_code,
                    'members': '\n'.join(members),
                    'members_init': ',\n'.join(members_init)
                    })
        else:
            raise ValueError('Only methods with one or two arguments allowed.')

        return code
