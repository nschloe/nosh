#ifndef NOSH_FVM_OPERATOR_H
#define NOSH_FVM_OPERATOR_H

#include "mesh.hpp"
#include "operator_core_boundary.hpp"
#include "operator_core_dirichlet.hpp"
#include "operator_core_edge.hpp"
#include "operator_core_vertex.hpp"

namespace nosh
{
  class fvm_operator:
    public Tpetra::Operator<double,int,int>
  {
    public:
      fvm_operator(
          const std::shared_ptr<const nosh::mesh> & _mesh,
          const std::set<std::shared_ptr<const operator_core_edge>> & operator_core_edges,
          const std::set<std::shared_ptr<const operator_core_vertex>> & operator_core_vertexs,
          const std::set<std::shared_ptr<const operator_core_boundary>> & operator_core_boundarys,
          const std::set<std::shared_ptr<const Tpetra::Operator<double,int,int>> & operators,
          const std::set<std::shared_ptr<const operator_core_dirichlet>> & dbcs
          ) :
        mesh(_mesh),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        apply_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: fvm_operator::apply")),
#endif
        operator_core_edges_(operator_core_edges),
        operator_core_vertexs_(operator_core_vertexs),
        operator_core_boundarys_(operator_core_boundarys),
        operators_(operators),
        dbcs_(dbcs)
        {
        }

      virtual
      ~fvm_operator()
      {};

      virtual
      apply(
          const MultiVector<double,int,int> & x,
          MultiVector<double,int,int> & y,
          Teuchos::ETransp mode=Teuchos::NO_TRANS,
          double alpha=Teuchos::ScalarTraits<double>::one(),
          dobule beta=Teuchos::ScalarTraits<double>::zero()
          ) const
      {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor tm(*apply_time_);
#endif
#ifndef NDEBUG
        TEUCHOS_ASSERT(this->mesh);
#endif
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            mode != Teuchos::NO_TRANS,
            "Only untransposed applies supported."
            );
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            alpha != 1.0,
            "Only alpha==1.0 supported."
            );
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            beta != 0.0,
            "Only beta==0.0 supported."
            );
        y.putScalar(0.0);

        const auto x_data = x.getDataNonConst();
        auto y_data = y.getDataNonConst();

        this->apply_edge_contributions_(x_data, y_data);
        this->apply_vertex_contributions_(x_data, y_data);
        this->apply_domain_boundary_contributions_(x_data, y_data);

        const auto yk = Tpetra::MultiVector<double,int,int>(y, Teuchos::Copy);
        for (const auto & operator: operators) {
          operator->apply(x, yk);
          y.update(1.0, yk, 1.0);
        }

        this->apply_dbcs_(x_data, y_data);

        return;
      }

    protected:
      void
      apply_edge_contributions_(
          const Teuchos::ArrayRCP<const double> & x_data,
          const Teuchos::ArrayRCP<double> & y_data
          )
      {
        for (const auto & core: this->operator_core_edges_) {
          for (const auto & subdomain_id: core->subdomain_ids) {
            const auto edge_data = this->mesh->get_edge_data();
            // this->meshset interior edges
            const auto edges = this->mesh->get_edges(subdomain_id);
            for (const auto edge: edges) {
              const auto verts = this->mesh->get_vertex_tuple(edge);
              const auto lid = this->mesh->local_index(edge);
              auto vals = core->eval(
                  this->mesh->get_coords(verts[0]),
                  this->mesh->get_coords(verts[1]),
                  edge_data[lid].length,
                  edge_data[lid].covolume,
                  // TODO lids
                  x_data[lid[0]],
                  x_data[lid[1]]
                  );

              y_data[lid[0]] += vals[0];
              y_data[lid[1]] += vals[1];
            }

            // this->meshset boundary edges
            const auto half_edges = this->mesh->get_edges(
                subdomain_id + "_halfedges"
                );
            for (const auto edge: half_edges) {
              const auto verts = this->mesh->get_vertex_tuple(edge);

              const auto lid = this->mesh->local_index(edge);
              auto vals = core->eval(
                  this->mesh->get_coords(verts[0]),
                  this->mesh->get_coords(verts[1]),
                  edge_data[lid].length,
                  edge_data[lid].covolume,
                  x_data[lid[0]],
                  x_data[lid[1]]
                  );

              // check which one of the two verts is in this->meshset
              if (this->mesh->contains(subdomain_id, {verts[0]})) {
                y_data[lid[0]] += vals[0];
              } else if (this->mesh->contains(subdomain_id, {verts[1]})) {
                y_data[lid[1]] += vals[1];
              } else {
                TEUCHOS_TEST_FOR_EXCEPT_MSG(
                    true,
                    "Neither of the two edge vertices is contained in the subdomain."
                    );
              }
            }
          }
        }
      }

      void
      apply_vertex_contributions(
          const Teuchos::ArrayRCP<const double> & x_data,
          const Teuchos::ArrayRCP<double> & y_data
          )
      {
        const auto & control_volumes = this->mesh->control_volumes();
        const auto c_data = control_volumes->getData();
        for (const auto & core: this->operator_core_vertexs_) {
          for (const auto & subdomain_id: core->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto & vertex: verts) {
              const auto k = this->mesh->local_index(vertex);
              y_data[k] += core->eval(
                  this->mesh->get_coords(vertex),
                  c_data[k],
                  x_data[k]
                  );
            }
          }
        }
      }

      void
      apply_domain_boundary_contributions(
          const Teuchos::ArrayRCP<const double> & x_data,
          const Teuchos::ArrayRCP<double> & y_data
          )
      {
        const auto surfs = this->mesh->boundary_surface_areas();
        for (const auto core: this->operator_core_boundarys_) {
          for (const auto & subdomain_id: core->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto vert: verts) {
              const auto k = this->mesh->local_index(vert);
              y_data[k] += core->eval(
                  this->mesh->get_coords(vert),
                  surfs[k],
                  x_data[i]
                  );
            }
          }
        }
      }

      void
      apply_dbcs(
          const Teuchos::ArrayRCP<const double> & x_data,
          const Teuchos::ArrayRCP<double> & y_data
          )
      {
        for (const auto & bc: this->dbcs_) {
          for (const auto & subdomain_id: bc->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto & vertex: verts) {
              const auto k = this->mesh->local_index(vert);
              y_data[k] = core->eval(
                this->mesh->get_coords(vertex),
                x_data[k]
                );
            }
          }
        }
      }

    public:
      const std::shared_ptr<const nosh::mesh> mesh;

    private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      const Teuchos::RCP<Teuchos::Time> apply_time_;
#endif
      const std::set<std::shared_ptr<const operator_core_edge>> operator_core_edges_;
      const std::set<std::shared_ptr<const operator_core_vertex>> operator_core_vertexs_;
      const std::set<std::shared_ptr<const operator_core_boundary>> operator_core_boundarys_;
      const std::set<std::shared_ptr<const Tpetra::Operator<double,int,int>>> operators_;
      const std::set<std::shared_ptr<const nosh::operator_core_dirichlet>> dbcs_;
  };
} // namespace nosh

#endif // NOSH_FVM_OPERATOR_H
