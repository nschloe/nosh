class ${name}:
  public nosh::matrix_core
{
  public:
    ${name}()${members_init} {};

    virtual ~${name}() {};

    virtual
      nosh::edge_data
      edge_contrib(
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

    virtual
      nosh::vertex_data
      vertex_contrib(
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
