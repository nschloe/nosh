class ${name}:
  public nosh::boundary_core
{
  public:
    ${name}()${members_init} {};

    virtual ~${name}() {};

    virtual
      nosh::boundary_data
      eval(
          const Eigen::Vector3d & x,
          const double surface_area
          ) const
      {
        ${db_body}
        return {
          ${db_coeff},
          ${db_affine}
          };
      }

  private:
    ${members_declare}
}; // class ${name}
