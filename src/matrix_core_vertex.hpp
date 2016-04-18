#ifndef NOSH_MATRIX_CORE_VERTEX_H
#define NOSH_MATRIX_CORE_VERTEX_H

#include <moab/Core.hpp>

namespace nosh
{
  struct vertex_data {
    double lhs;
    double rhs;
  };

  class matrix_core_vertex
  {
    public:
      explicit matrix_core_vertex(
          const std::set<std::string> & _subdomain_ids = {"everywhere"}
          ):
        subdomain_ids(_subdomain_ids)
        {};

      virtual ~matrix_core_vertex() {};

      virtual
      vertex_data
      eval(const moab::EntityHandle & vertex) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_MATRIX_CORE_VERTEX_H
