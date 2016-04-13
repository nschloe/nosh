#ifndef NOSH_OPERATOR_CORE_EDGE_H
#define NOSH_OPERATOR_CORE_EDGE_H

#include <Eigen/Dense>

#include <parameter_object.hpp>

namespace nosh
{
  class operator_core_edge: parameter_object
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
          const moab::EnityHandle & edge,
          const Teuchos::ArrayRCP<const double> & u
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_OPERATOR_CORE_EDGE_H
