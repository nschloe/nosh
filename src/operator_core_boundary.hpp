#ifndef NOSH_OPERATOR_CORE_BOUNDARY_H
#define NOSH_OPERATOR_CORE_BOUNDARY_H

#include <parameter_object.hpp>

namespace nosh
{
  class operator_core_boundary: public parameter_object
  {
    public:
      explicit operator_core_boundary(
          const std::set<std::string> & _subdomain_ids = {"boundary"}
          ):
        subdomain_ids(_subdomain_ids)
        {};

      virtual ~operator_core_boundary() {};

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

#endif // NOSH_OPERATOR_CORE_BOUNDARY_H
