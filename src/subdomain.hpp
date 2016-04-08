#ifndef NOSH_SUBDOMAIN_HPP
#define NOSH_SUBDOMAIN_HPP

namespace nosh {
  class subdomain
  {
    public:
      explicit subdomain(
          const std::string & _id,
          const bool _is_boundary_only
          ):
        id(_id),
        is_boundary_only(_is_boundary_only)
      {}

      virtual ~subdomain() {};

      virtual bool
      is_inside(const Eigen::Vector3d & x) const = 0;

    public:
      const std::string id;
      const bool is_boundary_only;
  };
}
#endif // NOSH_SUBDOMAIN_HPP
