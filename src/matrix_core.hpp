#ifndef NOSH_MATRIX_CORE_H
#define NOSH_MATRIX_CORE_H

namespace nosh
{
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
      std::vector<std::vector<double>>
      edge_contrib(
          const Eigen::Vector3d & x0,
          const Eigen::Vector3d & x1,
          const double edge_length,
          const double covolume
          ) const = 0;

      virtual
      double
      vertex_contrib(
          const Eigen::Vector3d & x,
          const double control_volume
          ) const = 0;
  };
} // namespace nosh

#endif // NOSH_MATRIX_CORE_H
