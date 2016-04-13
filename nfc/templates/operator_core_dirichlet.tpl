class ${name}:
  public nosh::operator_core_dirichlet
{
  public:
    ${name}(): ${init} {}

    virtual ~${name}() {}

    virtual double
    eval(
      const moab::EntityHandle & vertex,
      const Teuchos::ArrayRCP<const double> & u
    ) const {
      ${eval_body}return ${eval_return_value};
    }
}; // class ${name}
