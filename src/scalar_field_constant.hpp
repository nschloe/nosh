#ifndef NOSH_SCALARFIELD_CONSTANT_H_
#define NOSH_SCALARFIELD_CONSTANT_H_
// =============================================================================
#include <map>
#include <string>

#include <Teuchos_RCP.hpp>

#include "scalar_field_base.hpp"

#include "mesh.hpp"
// =============================================================================
namespace nosh
{
namespace scalar_field
{
class constant: public base
{
public:
  constant(
      const nosh::mesh & mesh,
      const double c,
      const std::string & param1_name = "",
      const double param1_init_value = 0.0
      );

  Tpetra::Vector<double,int,int>
  create_p_init_(const Tpetra::Map<int,int> & map);

  ~constant();

  //! Get the parameter names and intial values.
  virtual
  const std::map<std::string, double>
  get_parameters() const;

  virtual
  const Tpetra::Vector<double,int,int>
  get_v(const std::map<std::string, double> & params) const;

  virtual
  const Tpetra::Vector<double,int,int>
  get_dvdp(const std::map<std::string, double> & params,
          const std::string & param_name
        ) const;

protected:
private:
  const std::shared_ptr<const Tpetra::Map<int,int>> map_;
  const double c_;
  const std::string param1_name_;
  const double param1_init_value_;
};
} // namespace scalar_field
} // namespace nosh
#endif // NOSH_SCALARFIELD_CONSTANT_H_
