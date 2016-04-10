#ifndef NOSH_MATRIX_CORE_EDGE_H
#define NOSH_MATRIX_CORE_EDGE_H

#include <Eigen/Dense>

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
          const std::set<std::string> & _subdomain_ids = {"everywhere"}
          ):
        subdomain_ids(_subdomain_ids)
        {};

      virtual ~matrix_core_edge() {};

      virtual
      matrix_core_edge_data
      eval(
          const Eigen::Vector3d & x0,
          const Eigen::Vector3d & x1,
          const double edge_length,
          const double covolume
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_MATRIX_CORE_EDGE_H
