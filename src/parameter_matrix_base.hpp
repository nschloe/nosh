#ifndef NOSH_PARAMETER_MATRIX_VIRTUAL
#define NOSH_PARAMETER_MATRIX_VIRTUAL

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
namespace parameter_matrix
{
class base: public Tpetra::CrsMatrix<double,int,int>
{
public:
  explicit base(const std::shared_ptr<const nosh::mesh> &mesh);

  // Destructor.
  virtual
  ~base();

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

protected:
  const std::shared_ptr<const nosh::mesh> mesh_;

private:
  std::map<std::string, double> build_parameters_;
};
} // namespace parameter_matrix
} // namespace nosh

#endif // NOSH_PARAMETER_MATRIX_VIRTUAL
