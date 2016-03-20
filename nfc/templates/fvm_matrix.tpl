class ${name}:
  public nosh::fvm_matrix
{
  public:
    ${name}(
        ${constructor_args}
        ):
      ${members_init}
    {
      this->fill_();
    };

    virtual
      ~${name}()
      {};

  protected:
    virtual
      std::vector<std::vector<double>>
      edge_contrib(
          const Eigen::Vector3d & x0,
          const Eigen::Vector3d & x1,
          const double edge_length,
          const double edge_covolume
          ) const
      {
        ${edge_contrib_unused_args}
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
          const double control_volume
          ) const
      {
        ${vertex_contrib_unused_args}return ${vertex_contrib};
      }
  private:
    ${members_declare}
}; // class ${name}
