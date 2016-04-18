#include "scalar_field_explicit_values.hpp"

#include <map>
#include <string>

#include "mesh.hpp"

namespace nosh
{
namespace scalar_field
{
// ============================================================================
explicit_values::
explicit_values(const nosh::mesh & mesh,
               const std::string & field_name
             ) :
  node_values_(mesh.get_vector(field_name))
{
}
// ============================================================================
explicit_values::
~explicit_values()
{
}
// ============================================================================
const std::map<std::string, double>
explicit_values::
get_scalar_parameters() const
{
  std::map<std::string, double> m;
  m["beta"] = 1.0;
  return m;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
explicit_values::
get_v(const std::map<std::string, double> & params) const
{
  Tpetra::Vector<double,int,int> vals(*node_values_);
  // Scale by "beta"
  vals.scale(params.at("beta"));
  return vals;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
explicit_values::
get_dvdp(const std::map<std::string, double> & params,
        const std::string & param_name
      ) const
{
  (void) params;
  if (param_name.compare("beta") == 0) {
    return *node_values_;
  } else {
    return Tpetra::Vector<double,int,int>(
        node_values_->getMap(),
        true // zero out
        );
  }
}
// ============================================================================
} // namespace scalar_field
} // namespace nosh
