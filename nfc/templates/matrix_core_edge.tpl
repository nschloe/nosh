class ${name}:
  public nosh::matrix_core_edge
{
  public:
    ${name}()${members_init} {}

    virtual ~${name}() {}

    virtual
      nosh::matrix_core_edge_data
      eval(
          const Eigen::Vector3d & x0,
          const Eigen::Vector3d & x1,
          const double edge_length,
          const double edge_covolume
          ) const
      {
        ${edge_body}
        return {
          {
            {
              ${edge00},
                ${edge01},
            },
              {
                ${edge10},
                ${edge11},
              }
          },
          {
            ${edge_affine0},
            ${edge_affine1}
          }
        };
      }

  private:
    ${members_declare}
}; // class ${name}
