class ${name}:
  public nosh::fvm_matrix
{
  public:
    ${name}(
        const std::shared_ptr<const nosh::mesh> & _mesh
        ):
      nosh::fvm_matrix(
        _mesh,
        ${boundary_conditions}
        )
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
          const double edge_coefficient,
          const Eigen::Vector3d & edge_midpoint
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
}; // class ${name}
