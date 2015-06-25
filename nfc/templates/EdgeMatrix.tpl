class ${name}:
  public Nosh::EdgeMatrix
{
  public:
    ${name}(
        const std::shared_ptr<const Nosh::Mesh> & _mesh,
        const std::set<std::shared_ptr<const Nosh::DirichletBC>> & _bcs
        ):
      Nosh::EdgeMatrix(_mesh, _bcs)
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
          const double controlVolume0,
          const double controlVolume1
          ) const
      {
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
};
