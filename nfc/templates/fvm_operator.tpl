class ${name}:
  public nosh::fvm_operator
{
  public:
    ${name}(
        ${constructor_args}
        ):
      edge_cores_(${init_edge_cores}),
      vertex_cores_(${init_vertex_cores}),
      boundary_cores_(${init_boundary_cores}),
      dirichlets_(${init_dirichlets}),
      operators_(${init_operators}),
      nosh::fvm_operator(
        _mesh,
        edge_cores_,
        vertex_cores_,
        boundary_cores_,
        dirichlets_,
        operators_
        )
      ${members_init}
    {
    }

    virtual
      ~${name}()
      {}

    ${extra_methods}

  private:
    std::vector<std::shared_ptr<nosh::operator_core_edge>> edge_cores_;
    std::vector<std::shared_ptr<nosh::operator_core_vertex>> vertex_cores_;
    std::vector<std::shared_ptr<nosh::operator_core_boundary>> boundary_cores_;
    std::vector<std::shared_ptr<nosh::operator_core_dirichlet>> dirichlets_;
    std::vector<std::shared_ptr<Tpetra::Operator<double,int,int>>> operators_;
    ${members_declare}
}; // class ${name}
