#ifndef NOSH_PARAMETEROBJECT
#define NOSH_PARAMETEROBJECT

#include <map>
#include <string>

#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_RCP.hpp>

// forward declarations
namespace nosh
{
class mesh;
}

namespace nosh
{
class parameter_object
{
public:
  parameter_object();

  // Destructor.
  virtual
  ~parameter_object();

  //! Fill the matrix with the parameter entries as given in params.
  //! Includes some caching logic for params.
  virtual
  void
  set_parameters(const std::map<std::string, double> &params) final;

  //! Get parameter map with their initial values.
  virtual
  const std::map<std::string, double>
  get_parameters() const = 0;

protected:
  //! Fill the matrix with the parameter entries as given in params.
  virtual
  void
  refill_(const std::map<std::string, double> &params) = 0;

private:
  std::map<std::string, double> build_parameters_;
};
} // namespace nosh

#endif // NOSH_PARAMETEROBJECT
