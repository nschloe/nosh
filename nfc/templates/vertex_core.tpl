class ${name}:
  public nosh::vertex_core
{
  public:
    ${name}()${members_init} {}

    virtual ~${name}() {}

    virtual
      nosh::vertex_data
      eval(
          const Eigen::Vector3d & x,
          const double control_volume
          ) const
      {
        ${vertex_body}
        return {
          ${vertex_contrib},
          ${vertex_affine}
          };
      }

  private:
    ${members_declare}
}; // class ${name}
