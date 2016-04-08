#ifndef NOSH_VERTEXCORE_H
#define NOSH_VERTEXCORE_H

namespace nosh
{
  struct vertex_data {
    double lhs;
    double rhs;
  };

  class vertex_core
  {
    public:
      explicit vertex_core(
          const std::set<std::string> & _subdomain_ids = {"everywhere"}
          ):
        subdomain_ids(_subdomain_ids)
        {};

      virtual ~vertex_core() {};

      virtual
      vertex_data
      eval(
          const Eigen::Vector3d & x,
          const double control_volume
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_VERTEXCORE_H
