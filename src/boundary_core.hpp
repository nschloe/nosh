#ifndef NOSH_BOUNDARYCORE_H
#define NOSH_BOUNDARYCORE_H

namespace nosh
{
  struct boundary_data {
    double lhs;
    double rhs;
  };

  class boundary_core
  {
    public:
      explicit boundary_core(const std::set<std::string> & _subdomain_ids):
        subdomain_ids(_subdomain_ids)
        {};

      virtual ~boundary_core() {};

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

#endif // NOSH_BOUNDARYCORE_H
