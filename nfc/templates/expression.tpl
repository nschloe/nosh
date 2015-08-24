class ${name}:
  public nosh::expression
{
  public:
    ${name}():
      nosh::expression(${degree})
    {}

    virtual
    ~${name}()
    {}

    virtual
    double
    eval(const Eigen::Vector3d & x) const
    {
      return ${eval};
    };
};
