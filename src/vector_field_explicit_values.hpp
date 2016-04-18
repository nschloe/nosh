#ifndef NOSH_VECTORFIELD_EXPLICITVALUES_H_
#define NOSH_VECTORFIELD_EXPLICITVALUES_H_

#include <map>
#include <string>

#include "vector_field_base.hpp"
#include "mesh.hpp"

namespace nosh
{
namespace vector_field
{
class explicit_values : public base
{
public:
  explicit_values(
      const nosh::mesh &mesh,
      const std::string &field_name,
      const double mu
      );

  virtual
  ~explicit_values();

  virtual
  void
  set_parameters(const std::map<std::string, double> & params);

  //! Get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  get_scalar_parameters() const;

  double
  get_edge_projection(const unsigned int edge_index) const;

  double
  get_d_edge_projection_dp(
      const unsigned int edge_index,
      const std::string & dParamName
      ) const;

  virtual
  Eigen::Vector3d
  eval(const Eigen::Vector3d & x) const {
    // mu * 0.5  * X.cross(B)
    return mu_ * 0.5 * Eigen::Vector3d({-x[1], x[0], 0.0});
  };

  virtual
  unsigned int
  degree() const {
    return 1;
  };

protected:
private:
  double mu_;

  std::vector<double> edgeProjectionCache_;
};
} // namespace vector_field
} // namespace nosh
#endif // NOSH_VECTORFIELD_EXPLICITVALUES_H_
