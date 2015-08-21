class ${name}:
  public Nosh::DirichletBC
{
  public:
    virtual bool
    isInside(const Eigen::Vector3d & x) const {
      ${inside_void}return ${insideCondition};
    }
    virtual double
    eval(const Eigen::Vector3d & x) const {
      ${eval_void}return ${evalReturnValue};
    }
};
