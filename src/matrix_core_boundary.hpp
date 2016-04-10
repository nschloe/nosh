#ifndef NOSH_MATRIX_CORE_BOUNDARY_H
#define NOSH_MATRIX_CORE_BOUNDARY_H

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
          const std::set<std::string> & _subdomain_ids = {"boundary"}
          ):
        subdomain_ids(_subdomain_ids)
        {};

      virtual ~matrix_core_boundary() {};

      virtual
      boundary_data
      eval(
          const Eigen::Vector3d & x,
          const double surface_area
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_MATRIX_CORE_BOUNDARY_H
