#include "scalar_field_constant.hpp"

#include <map>
#include <string>

#include <Tpetra_Map.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

namespace nosh
{
namespace scalar_field
{
// ============================================================================
constant::
constant(
    const nosh::mesh & mesh,
    const double c,
    std::string  param1_name,
    const double param1_init_value
    ):
  map_(mesh.map()),
  c_(c),
  param1_name_(std::move(param1_name)),
  param1_init_value_(param1_init_value)
{
}
// ============================================================================
constant::
~constant()
= default;
// ============================================================================
const std::map<std::string, double>
constant::
get_scalar_parameters() const
{
  std::map<std::string, double> m;
  if (!param1_name_.empty()) {
    m[param1_name_] = param1_init_value_;
}
  return m;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
constant::
get_v(const std::map<std::string, double> & params) const
{
  // Create constant-valued vector.
  Tpetra::Vector<double,int,int> vals(Teuchos::rcp(map_));

  auto it = params.find(param1_name_);
  if (it != params.end()) {
    vals.putScalar(c_ + it->second);
  } else {
    vals.putScalar(c_);
  }

  return vals;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
constant::
get_dvdp(
    const std::map<std::string, double> & params,
    const std::string & param_name
    ) const
{
  (void) params;
  // Create zeroed-out vector.
  Tpetra::Vector<double,int,int> vals(Teuchos::rcp(map_), true);
  if (param_name.compare(param1_name_) == 0) {
    vals.putScalar(1.0);
  }

  return vals;
}
// ============================================================================
} // namespace scalar_field
} // namespace nosh
