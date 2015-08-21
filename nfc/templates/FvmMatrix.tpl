class ${name}:
  public Nosh::FvmMatrix
{
  public:
    ${name}(
        const std::shared_ptr<const Nosh::Mesh> & _mesh,
        const std::set<std::shared_ptr<const Nosh::DirichletBC>> & _bcs
        ):
      Nosh::FvmMatrix(_mesh, _bcs)
    {
      this->fill_();
    };

    virtual
      ~${name}()
      {};

  protected:
    virtual
      std::vector<std::vector<double>>
      edgeContrib(
          const double edgeCoefficient,
          const Eigen::Vector3d & edgeMidpoint
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
      vertexContrib(
          const double controlVolume
          ) const
      {
        ${vertex_contrib_unused_args}
        return ${vertex_contrib};
      }
  private:
};
