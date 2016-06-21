#ifndef NOSH_MATRIX_CORE_EDGE_H
#define NOSH_MATRIX_CORE_EDGE_H

#include <set>
#include <Eigen/Dense>
#include <moab/Core.hpp>

namespace nosh
{
  struct matrix_core_edge_data {
    std::vector<std::vector<double>> lhs;
    std::vector<double> rhs;
  };

  class matrix_core_edge
  {
    public:
      explicit matrix_core_edge(
          std::set<std::string>  _subdomain_ids = {"everywhere"}
          ):
        subdomain_ids(std::move(_subdomain_ids))
        {};

      virtual ~matrix_core_edge() = default;

      virtual
      matrix_core_edge_data
      eval(const moab::EntityHandle & edge) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_MATRIX_CORE_EDGE_H
