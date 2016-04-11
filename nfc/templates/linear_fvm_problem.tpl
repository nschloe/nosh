class ${name}:
  public nosh::linear_problem
{
  public:
    ${name}(
        ${constructor_args}
        ):
      ${members_init}
    {
      this->fill();
    }

    virtual
      ~${name}()
      {}

  private:
    ${members_declare}
}; // class ${name}
