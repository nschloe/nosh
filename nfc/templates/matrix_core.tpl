class ${name}:
  public nosh::matrix_core
{
  public:
    ${name}()${members_init} {};

    virtual ~${name}() {};

    virtual
      std::vector<std::vector<double>>
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
            ${edge00},
            ${edge01},
          },
          {
            ${edge10},
            ${edge11},
          }
        };
      }

    virtual
      double
      vertex_contrib(
          const Eigen::Vector3d & x,
          const double control_volume
          ) const
      {
        ${vertex_body}
        return ${vertex_contrib};
      }

  private:
    ${members_declare}
}; // class ${name}
