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

    ${extra_methods}

  private:
    ${members_declare}
}; // class ${name}
