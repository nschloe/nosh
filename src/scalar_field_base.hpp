#ifndef NOSH_SCALARFIELD_BASE_H_
#define NOSH_SCALARFIELD_BASE_H_
// =============================================================================
#include <map>
#include <string>

#include <Tpetra_Vector.hpp>

namespace nosh
{
namespace scalar_field
{
class base
{
public:
  virtual
  ~base() = default;

  virtual
  const Tpetra::Vector<double,int,int>
  get_v(const std::map<std::string, double> & params) const = 0;

  virtual
  const Tpetra::Vector<double,int,int>
  get_dvdp(
      const std::map<std::string, double> & params,
      const std::string & param_name
      ) const = 0;

  //! Get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  get_scalar_parameters() const = 0;

protected:
private:
};
} // namespace scalar_field
} // namespace nosh
#endif // NOSH_SCALARFIELD_BASE_H_
