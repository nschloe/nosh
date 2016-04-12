#ifndef NOSH_PARAMETEROBJECT
#define NOSH_PARAMETEROBJECT

#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_RCP.hpp>

#include <map>
#include <string>

namespace nosh
{
class parameter_object
{
public:
  parameter_object():
    build_parameters_scalar_(),
    build_parameters_vector_()
  {
  }

  virtual
  ~parameter_object() {};

  //! Fill the matrix with the parameter entries as given in params.
  //! Includes some caching logic for params.
  virtual
  void
  set_parameters(
      const std::map<std::string, double> & scalar_params,
      const std::map<std::string, std::shared_ptr<Tpetra::Vector<double, int, int>>> & vector_params
      ) final;

  //! Get scalar parameter map with their initial values.
  virtual
  const std::map<std::string, double>
  get_scalar_parameters() const
  {
    return {};
  };

  //! Get vector parameter map with their initial values.
  virtual
  std::map<std::string, std::shared_ptr<Tpetra::Vector<double, int, int>>>
  get_vector_parameters() const
  {
    return {};
  };

protected:
  //! Fill the matrix with the parameter entries as given in params.
  virtual
  void
  refill_(
      const std::map<std::string, double> & scalar_params,
      const std::map<std::string, std::shared_ptr<Tpetra::Vector<double, int, int>>> & vector_params
      )
  {
    // By default, don't do anything
    (void) scalar_params;
    (void) vector_params;
  }

private:
  std::map<std::string, double> build_parameters_scalar_;
  std::map<std::string, std::shared_ptr<Tpetra::Vector<double, int, int>>> build_parameters_vector_;
};
}  // namespace nosh
#endif  // NOSH_PARAMETEROBJECT
