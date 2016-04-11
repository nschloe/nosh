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
        subdomain_ids(_subdomain_ids),
        scalar_parameters_({}),
        vector_parameters_({})
        {};

      virtual ~operator_core_boundary() {};

      virtual
      double
      eval(
          const Eigen::Vector3d & x,
          const double surface_area,
          const double u
          ) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
      const std::map<std::string, double> scalar_parameters_;
      const std::map<std::string, std::shared_ptr<Tpetra::Vector<double,int,int>>> vector_parameters_;
  };
} // namespace nosh

#endif // NOSH_OPERATOR_CORE_BOUNDARY_H
