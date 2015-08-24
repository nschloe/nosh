class ${name}:
  public nosh::dirichlet_bc
{
  public:
    virtual bool
    is_inside(const Eigen::Vector3d & x) const {
      ${inside_void}return ${inside_condition};
    }
    virtual double
    eval(const Eigen::Vector3d & x) const {
      ${eval_void}return ${eval_return_value};
    }
};
