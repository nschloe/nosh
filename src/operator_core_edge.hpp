#ifndef NOSH_OPERATOR_CORE_EDGE_H
#define NOSH_OPERATOR_CORE_EDGE_H

#include <Eigen/Dense>

namespace nosh
{
  class operator_core_edge
  {
    public:
      explicit operator_core_edge(
          const std::set<std::string> & _subdomain_ids = {"everywhere"}
          ):
        subdomain_ids(_subdomain_ids)
        {};

      virtual ~operator_core_edge() {};

      virtual
      std::tuple<double,double>
      eval(
          const Eigen::Vector3d & x0,
          const Eigen::Vector3d & x1,
          const double edge_length,
          const double covolume,
          const double u0,
          const double u1
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_OPERATOR_CORE_EDGE_H
