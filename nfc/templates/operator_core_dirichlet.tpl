class ${name}:
  public nosh::operator_core_dirichlet
{
  public:
    ${name}(): ${init} {}

    virtual ~${name}() {}

    virtual double
    eval(
      const Eigen::Vector3d & x,
      const double u
    ) const {
      ${eval_body}return ${eval_return_value};
    }
}; // class ${name}
