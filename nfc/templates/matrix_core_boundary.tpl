class ${name}: public nosh::matrix_core_boundary
{
  public:
    ${name}()${members_init} {}

    virtual ~${name}() {}

    virtual
      nosh::boundary_data
      eval(
          const Eigen::Vector3d & x,
          const double surface_area
          ) const
      {
        ${body}
        return {
          ${coeff},
          ${affine}
          };
      }

  private:
    ${members_declare}
}; // class ${name}
