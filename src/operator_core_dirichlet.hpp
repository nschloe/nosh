#ifndef NOSH_OPERATOR_CORE_DIRICHLET_HPP
#define NOSH_OPERATOR_CORE_DIRICHLET_HPP

#include <Eigen/Dense>
#include <moab/Core.hpp>

#include "parameter_object.hpp"

namespace nosh {
  class operator_core_dirichlet: public parameter_object
  {
    public:
      explicit operator_core_dirichlet(
          std::set<std::string>  _subdomain_ids
          ):
        subdomain_ids(std::move(_subdomain_ids))
      {};

      ~operator_core_dirichlet() override = default;

      virtual
      double
      eval(
          const moab::EntityHandle & vertex,
          const Teuchos::ArrayRCP<const double> & u
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
}  // namespace nosh
#endif // NOSH_DIRICHLET_BC_HPP
