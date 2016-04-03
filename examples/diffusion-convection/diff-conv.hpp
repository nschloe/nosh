#ifndef POISSON_HPP
#define POISSON_HPP

#include <nosh.hpp>

namespace diff_conv
{
class bc1:
  public nosh::dirichlet_bc
{
public:
  virtual bool
  is_inside(const Eigen::Vector3d & x) const
  {
    (void) x;
    return true;
  }
  virtual double
  eval(const Eigen::Vector3d & x) const
  {
    (void) x;
    return 0.0;
  }
}; // class bc1

class core0:
  public nosh::matrix_core
{
public:
  core0() {};

  virtual ~core0() {};

  virtual
  std::vector<std::vector<double>>
                                edge_contrib(
                                  const Eigen::Vector3d & x0,
                                  const Eigen::Vector3d & x1,
                                  const double edge_length,
                                  const double edge_covolume
                                ) const
  {
    (void) x0;
    (void) x1;
    const double alpha = edge_covolume/edge_length;
    const auto n = (x1 - x0) / edge_length;
    const Eigen::Vector3d a(1, 1, 0);
    const double beta = edge_covolume / 2 * n.dot(a);
    return
    {
      {
        alpha + beta,
        -alpha + beta,
      },
      {
        -alpha - beta,
        alpha - beta,
      }
    };
  }

  virtual
  double
  vertex_contrib(
    const Eigen::Vector3d & x,
    const double control_volume
  ) const
  {
    (void) control_volume;
    (void) x;
    return 0.0;
  }

private:

}; // class core0

class laplace:
  public nosh::fvm_matrix
{
public:
  laplace(
    const std::shared_ptr<const nosh::mesh> & _mesh
  ):
    nosh::fvm_matrix(
        _mesh,
        {std::make_shared<core0>()},
        {std::make_shared<bc1>()}
        )
  {
    this->fill_();
  };

  virtual
  ~laplace()
  {};

private:

}; // class laplace

class f:
  public nosh::expression
{
public:
  f():
    nosh::expression(0)
  {}

  virtual
  ~f()
  {}

  virtual
  double
  operator()(const Eigen::Vector3d & x) const
  {
    (void) x;
    return 1.0;
  };
};

} // namespace diff_conv

#endif // POISSON_HPP
