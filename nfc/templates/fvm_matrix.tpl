class ${name}:
  public nosh::fvm_matrix
{
  public:
    ${name}(
        const std::shared_ptr<const nosh::mesh> & _mesh,
        const std::set<std::shared_ptr<const nosh::dirichlet_bc>> & _bcs
        ):
      nosh::fvm_matrix(_mesh, _bcs)
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
