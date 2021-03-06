#ifndef NOSH_MATRIX_CORE_BOUNDARY_H
#define NOSH_MATRIX_CORE_BOUNDARY_H

#include <set>
#include <moab/Core.hpp>

namespace nosh
{
  struct boundary_data {
    double lhs;
    double rhs;
  };

  class matrix_core_boundary
  {
    public:
      explicit matrix_core_boundary(
          std::set<std::string>  _subdomain_ids = {"boundary"}
          ):
        subdomain_ids(std::move(_subdomain_ids))
        {};

      virtual ~matrix_core_boundary() = default;

      virtual
      boundary_data
      eval(const moab::EntityHandle & vertex) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_MATRIX_CORE_BOUNDARY_H
