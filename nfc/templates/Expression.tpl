class ${name}:
  public Nosh::Expression
{
  public:
    ${name}():
      Nosh::Expression(${degree})
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
