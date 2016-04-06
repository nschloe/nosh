#ifndef NOSH_MATRIX_CORE_H
#define NOSH_MATRIX_CORE_H

namespace nosh
{
  struct edge_data {
    std::vector<std::vector<double>> lhs;
    std::vector<double> rhs;
  };

  struct vertex_data {
    double lhs;
    double rhs;
  };

  class matrix_core
  {
    public:
      matrix_core() {};
      virtual ~matrix_core() {};

      virtual
      bool
      is_inside(
          const Eigen::Vector3d & x
          ) const
      {
        (void) x;
        return true;
      };

      virtual
      edge_data
      edge_contrib(
          const Eigen::Vector3d & x0,
          const Eigen::Vector3d & x1,
          const double edge_length,
          const double covolume
          ) const = 0;

      virtual
      vertex_data
      vertex_contrib(
          const Eigen::Vector3d & x,
          const double control_volume
          ) const = 0;
  };
} // namespace nosh

#endif // NOSH_MATRIX_CORE_H
