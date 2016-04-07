#ifndef NOSH_EDGECORE_H
#define NOSH_EDGECORE_H

namespace nosh
{
  struct edge_core_data {
    std::vector<std::vector<double>> lhs;
    std::vector<double> rhs;
  };

  class edge_core
  {
    public:
      explicit edge_core(const std::set<std::string> & _subdomain_ids):
        subdomain_ids(_subdomain_ids)
        {};

      virtual ~edge_core() {};

      virtual
      edge_core_data
      eval(
          const Eigen::Vector3d & x0,
          const Eigen::Vector3d & x1,
          const double edge_length,
          const double covolume
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
} // namespace nosh

#endif // NOSH_EDGECORE_H
