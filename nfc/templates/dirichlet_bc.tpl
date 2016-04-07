class ${name}:
  public nosh::dirichlet_bc
{
  public:
    ${name}(): ${init} {}

    virtual ~${name}() {}

    virtual double
    eval(const Eigen::Vector3d & x) const {
      ${eval_void}return ${eval_return_value};
    }
}; // class ${name}
