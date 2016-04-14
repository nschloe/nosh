#ifndef NOSH_FVM_OPERATOR_H
#define NOSH_FVM_OPERATOR_H

#include <tuple>

#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

#include "mesh.hpp"
#include "operator_core_boundary.hpp"
#include "operator_core_dirichlet.hpp"
#include "operator_core_edge.hpp"
#include "operator_core_vertex.hpp"
#include "parameter_object.hpp"

namespace nosh
{
  class fvm_operator:
    public Tpetra::Operator<double,int,int>, public parameter_object
  {
    public:
      fvm_operator(
        const std::shared_ptr<const nosh::mesh> & _mesh,
        const std::vector<std::shared_ptr<operator_core_edge>> & edge_cores,
        const std::vector<std::shared_ptr<operator_core_vertex>> & vertex_cores,
        const std::vector<std::shared_ptr<operator_core_boundary>> & boundary_cores,
        const std::vector<std::shared_ptr<operator_core_dirichlet>> & dirichlets,
        const std::vector<std::shared_ptr<Tpetra::Operator<double,int,int>>> & operators
        ) :
      mesh(_mesh),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      apply_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: fvm_operator::apply")),
#endif
      edge_cores_(edge_cores),
      vertex_cores_(vertex_cores),
      boundary_cores_(boundary_cores),
      dirichlets_(dirichlets),
      operators_(operators)
      {
      }

      virtual
      ~fvm_operator()
      {};

      virtual
      void
      apply(
          const Tpetra::MultiVector<double,int,int> & x,
          Tpetra::MultiVector<double,int,int> & y,
          Teuchos::ETransp mode = Teuchos::NO_TRANS,
          double alpha = Teuchos::ScalarTraits<double>::one(),
          double beta = Teuchos::ScalarTraits<double>::zero()
          ) const
      {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor tm(*apply_time_);
#endif
#ifndef NDEBUG
        TEUCHOS_ASSERT(this->mesh);
#endif
        TEUCHOS_ASSERT_EQUALITY(x.getNumVectors(), 0);
        TEUCHOS_ASSERT_EQUALITY(y.getNumVectors(), 0);

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

        const auto x_data = x.getData(0);
        auto y_data = y.getDataNonConst(0);

        this->apply_edge_contributions_(x_data, y_data);
        this->apply_vertex_contributions_(x_data, y_data);
        this->apply_domain_boundary_contributions_(x_data, y_data);

        auto yk = Tpetra::MultiVector<double,int,int>(y, Teuchos::Copy);
        for (const auto & op: this->operators_) {
          op->apply(x, yk);
          y.update(1.0, yk, 1.0);
        }

        this->apply_dirichlets_(x_data, y_data);

        return;
      }

      virtual
      Teuchos::RCP<const Tpetra::Map<int,int>>
      getDomainMap() const
      {
        return Teuchos::rcp(this->mesh->map());
      }

      virtual
      Teuchos::RCP<const Tpetra::Map<int,int>>
      getRangeMap() const
      {
        return Teuchos::rcp(this->mesh->map());
      }


    protected:
      void
      apply_edge_contributions_(
          const Teuchos::ArrayRCP<const double> & x_data,
          const Teuchos::ArrayRCP<double> & y_data
          ) const
      {
        for (const auto & core: this->edge_cores_) {
          for (const auto & subdomain_id: core->subdomain_ids) {
            const auto edge_data = this->mesh->get_edge_data();
            // this->meshset interior edges
            const auto edges = this->mesh->get_edges(subdomain_id);
            for (const auto edge: edges) {
              auto vals = core->eval(edge, x_data);

              const auto verts = this->mesh->get_vertex_tuple(edge);
              const auto i0 = this->mesh->local_index(verts[0]);
              const auto i1 = this->mesh->local_index(verts[1]);
              y_data[i0] += std::get<0>(vals);
              y_data[i1] += std::get<1>(vals);
            }

            // this->meshset boundary edges
            const auto half_edges = this->mesh->get_edges(
                subdomain_id + "_halfedges"
                );
            for (const auto edge: half_edges) {
              auto vals = core->eval(edge, x_data);

              const auto verts = this->mesh->get_vertex_tuple(edge);
              // check which one of the two verts is in this->meshset
              if (this->mesh->contains(subdomain_id, {verts[0]})) {
                const auto i0 = this->mesh->local_index(verts[0]);
                y_data[i0] += std::get<0>(vals);
              } else if (this->mesh->contains(subdomain_id, {verts[1]})) {
                const auto i1 = this->mesh->local_index(verts[1]);
                y_data[i1] += std::get<1>(vals);
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
      apply_vertex_contributions_(
          const Teuchos::ArrayRCP<const double> & x_data,
          const Teuchos::ArrayRCP<double> & y_data
          ) const
      {
        const auto & control_volumes = this->mesh->control_volumes();
        const auto c_data = control_volumes->getData();
        for (const auto & core: this->vertex_cores_) {
          for (const auto & subdomain_id: core->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto & vertex: verts) {
              const auto k = this->mesh->local_index(vertex);
              y_data[k] += core->eval(vertex, x_data);
            }
          }
        }
      }

      void
      apply_domain_boundary_contributions_(
          const Teuchos::ArrayRCP<const double> & x_data,
          const Teuchos::ArrayRCP<double> & y_data
          ) const
      {
        const auto surfs = this->mesh->boundary_surface_areas();
        for (const auto core: this->boundary_cores_) {
          for (const auto & subdomain_id: core->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto vert: verts) {
              const auto k = this->mesh->local_index(vert);
              y_data[k] += core->eval(vert, x_data);
            }
          }
        }
      }

      void
      apply_dirichlets_(
          const Teuchos::ArrayRCP<const double> & x_data,
          const Teuchos::ArrayRCP<double> & y_data
          ) const
      {
        for (const auto & bc: this->dirichlets_) {
          for (const auto & subdomain_id: bc->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto & vertex: verts) {
              const auto k = this->mesh->local_index(vertex);
              y_data[k] = bc->eval(vertex, x_data);
            }
          }
        }
      }

    public:
      const std::shared_ptr<const nosh::mesh> mesh;

    protected:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      const Teuchos::RCP<Teuchos::Time> apply_time_;
#endif
      const std::vector<std::shared_ptr<operator_core_edge>> edge_cores_;
      const std::vector<std::shared_ptr<operator_core_vertex>> vertex_cores_;
      const std::vector<std::shared_ptr<operator_core_boundary>> boundary_cores_;
      const std::vector<std::shared_ptr<operator_core_dirichlet>> dirichlets_;
      const std::vector<std::shared_ptr<Tpetra::Operator<double,int,int>>> operators_;
  };
} // namespace nosh

#endif // NOSH_FVM_OPERATOR_H
