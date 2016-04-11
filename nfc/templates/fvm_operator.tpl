class ${name}:
  public nosh::fvm_operator
{
  public:
    ${name}(
        ${constructor_args}
        ):
      ${members_init}
    {
    }

    virtual
      ~${name}()
      {}

  private:
    ${members_declare}
}; // class ${name}
