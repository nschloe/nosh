class ${name}:
  public Nosh::DirichletBC
{
  public:
    virtual bool
    isInside(const Eigen::Vector3d & x) const {
      return ${insideCondition};
    }
    virtual double
    eval(const Eigen::Vector3d & x) const {
      return ${evalReturnValue};
    }
};
