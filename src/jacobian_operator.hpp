#ifndef NOSH_JACOBIANOPERATOR_H
#define NOSH_JACOBIANOPERATOR_H

#include <map>
#include <string>

#include <Tpetra_Vector.hpp>
#include <Tpetra_Operator.hpp>
#include <Teuchos_RCP.hpp>

#include "parameter_matrix_keo.hpp"

// forward declarations
namespace nosh
{
  class mesh;
  namespace scalar_field
  {
    class base;
  }
} // namespace nosh

namespace nosh
{

class jacobian_operator : public Tpetra::Operator<double,int,int>
{
public:
  jacobian_operator(
      const std::shared_ptr<const nosh::mesh> &mesh,
      const std::shared_ptr<const nosh::scalar_field::base> &scalar_potential,
      const std::shared_ptr<const nosh::scalar_field::base> &thickness,
      const std::shared_ptr<nosh::parameter_matrix::keo> &keo
      );

  // Destructor.
  ~jacobian_operator ();

  virtual void
  apply(
      const Tpetra::MultiVector<double,int,int> &X,
      Tpetra::MultiVector<double,int,int> &Y,
      Teuchos::ETransp mode,
      double alpha,
      double beta
      ) const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>> getDomainMap() const
  {
    return keo_->getDomainMap();
  }

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>> getRangeMap() const
  {
    return keo_->getRangeMap();
  }

public:
  void
  rebuild(
      const std::map<std::string, double> & params,
      const Tpetra::Vector<double,int,int> & current_x
      );

protected:
private:
  void
  rebuild_diags_(
      const std::map<std::string, double> & params,
      const Tpetra::Vector<double,int,int> & current_x
      );

private:
  const std::shared_ptr<const nosh::mesh> mesh_;
  const std::shared_ptr<const nosh::scalar_field::base> scalar_potential_;
  const std::shared_ptr<const nosh::scalar_field::base> thickness_;

  const std::shared_ptr<nosh::parameter_matrix::keo> keo_;
  Tpetra::Vector<double,int,int> diag0_;
  Tpetra::Vector<double,int,int> diag1b_;
};
} // namespace nosh

#endif // NOSH_JACOBIANOPERATOR_H
