#ifndef NOSH_MATRIX_CORE_DIRICHLET_HPP
#define NOSH_MATRIX_CORE_DIRICHLET_HPP

#include <set>
#include <Eigen/Dense>
#include <moab/Core.hpp>

namespace nosh {
  class matrix_core_dirichlet
  {
    public:
      explicit matrix_core_dirichlet(
          std::set<std::string>  _subdomain_ids
          ):
        subdomain_ids(std::move(_subdomain_ids))
      {};

      virtual ~matrix_core_dirichlet() = default;

      virtual
      double
      eval(const moab::EntityHandle & vertex) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
}  // namespace nosh
#endif // NOSH_DIRICHLET_BC_HPP
