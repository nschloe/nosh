#ifndef NOSH_SCALARFIELD_EXPLICITVALUES_H_
#define NOSH_SCALARFIELD_EXPLICITVALUES_H_

#include <map>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Tpetra_Vector.hpp>

#include "scalar_field_base.hpp"
#include "mesh.hpp"

namespace nosh
{
namespace scalar_field
{
class explicit_values : public base
{
public:
  explicit_values(
      const nosh::mesh &mesh,
      const std::string &field_name
      );

  virtual
  ~explicit_values();

//! Get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  get_parameters() const;

  virtual
  const Tpetra::Vector<double,int,int>
  get_v(const std::map<std::string, double> & params) const;

  virtual
  const Tpetra::Vector<double,int,int>
  get_dvdp(
      const std::map<std::string, double> & params,
      const std::string & param_name
      ) const;

protected:
private:
  const std::shared_ptr<const Tpetra::Vector<double,int,int>> node_values_;
};
} // namespace scalar_field
} // namespace nosh
#endif // NOSH_SCALARFIELD_EXPLICITVALUES_H_
