class ${name}:
  public nosh::matrix_core_dirichlet
{
  public:
    ${name}(): ${init} {}

    virtual ~${name}() {}

    virtual double
    eval(const Eigen::Vector3d & x) const {
      ${eval_body}return ${eval_return_value};
    }
}; // class ${name}
