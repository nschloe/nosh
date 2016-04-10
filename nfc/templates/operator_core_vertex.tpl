class ${name}:
  public nosh::operator_core_vertex
{
  public:
    ${name}()${members_init} {}

    virtual ~${name}() {}

    virtual
      double
      eval(
          const Eigen::Vector3d & x,
          const double control_volume,
          const double u0
          ) const
      {
        ${vertex_body}
        return ${return_value};
      }

  private:
    ${members_declare}
}; // class ${name}
