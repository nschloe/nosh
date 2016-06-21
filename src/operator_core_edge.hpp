#ifndef NOSH_OPERATOR_CORE_EDGE_H
#define NOSH_OPERATOR_CORE_EDGE_H

#include <Eigen/Dense>
#include <moab/Core.hpp>

#include "parameter_object.hpp"

namespace nosh
{
  class operator_core_edge: parameter_object
  {
    public:
      explicit operator_core_edge(
          std::set<std::string>  _subdomain_ids = {"everywhere"}
          ):
        subdomain_ids(std::move(_subdomain_ids))
        {};

      ~operator_core_edge() override = default;

      virtual
      std::tuple<double,double>
      eval(
          const moab::EntityHandle & edge,
          const Teuchos::ArrayRCP<const double> & u
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_OPERATOR_CORE_EDGE_H
