#ifndef NOSH_VECTORFIELD_BASE_H_
#define NOSH_VECTORFIELD_BASE_H_

#include <string>
#include <map>

#include <Eigen/Dense>

namespace nosh
{
namespace vector_field
{
class base
{
public:
  virtual
  ~base(){};

  virtual
  void
  set_parameters(const std::map<std::string, double> & params) = 0;

  //! get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  get_scalar_parameters() const = 0;

  //! Projection of the vector field onto an edge at the midpoint of the edge.
  virtual
  double
  get_edge_projection(const unsigned int edge_index) const = 0;

  virtual
  double
  get_d_edge_projection_dp(
      const unsigned int edge_index,
      const std::string & dParamName
      ) const = 0;

  virtual
  Eigen::Vector3d
  eval(const Eigen::Vector3d & x) const = 0;

  virtual
  unsigned int
  degree() const = 0;

protected:
private:
};
} // namespace vector_field
} // namespace nosh
#endif // NOSH_VECTORFIELD_BASE_H_
