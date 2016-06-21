#ifndef NOSH_OPERATOR_CORE_VERTEX_H
#define NOSH_OPERATOR_CORE_VERTEX_H

#include "parameter_object.hpp"

#include <moab/Core.hpp>

namespace nosh
{
  class operator_core_vertex: public parameter_object
  {
    public:
      explicit operator_core_vertex(
          std::set<std::string>  _subdomain_ids = {"everywhere"}
          ):
        subdomain_ids(std::move(_subdomain_ids))
        {};

      ~operator_core_vertex() override = default;

      virtual
      double
      eval(
          const moab::EntityHandle & vertex,
          const Teuchos::ArrayRCP<const double> & u
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_OPERATOR_CORE_VERTEX_H
