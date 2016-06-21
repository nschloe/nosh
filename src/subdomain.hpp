#ifndef NOSH_SUBDOMAIN_HPP
#define NOSH_SUBDOMAIN_HPP

#include <Eigen/Dense>

namespace nosh {
  class subdomain
  {
    public:
      explicit subdomain(
          std::string _id,
          const bool _is_boundary_only
          ):
        id(std::move(_id)),
        is_boundary_only(_is_boundary_only)
      {}

      virtual ~subdomain() = default;

      virtual bool
      is_inside(const Eigen::Vector3d & x) const = 0;

    public:
      const std::string id;
      const bool is_boundary_only;
  };
}  // namespace nosh
#endif // NOSH_SUBDOMAIN_HPP
