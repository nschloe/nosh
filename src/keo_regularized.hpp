#ifndef NOSH_KEOPRECONDITIONER_H
#define NOSH_KEOPRECONDITIONER_H

#include <map>
#include <string>

#include <Tpetra_Vector.hpp>
#include <Tpetra_Operator.hpp>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Teuchos_FancyOStream.hpp>

#include <MueLu_TpetraOperator.hpp>

// forward declarations
namespace nosh
{
  class mesh;
  namespace scalar_field
  {
    class base;
  }
  namespace vector_field
  {
    class base;
  }
  namespace parameter_matrix
  {
    class keo;
  }
} // namespace nosh

namespace nosh
{
class keo_regularized : public Tpetra::Operator<double,int,int>
{
public:
  keo_regularized(
      const std::shared_ptr<const nosh::mesh> &mesh,
      const std::shared_ptr<const nosh::scalar_field::base> &thickness,
      const std::shared_ptr<nosh::vector_field::base> &mvp
      );

  // Destructor.
  ~keo_regularized();

  virtual void
  apply(
      const Tpetra::MultiVector<double,int,int> &X,
      Tpetra::MultiVector<double,int,int> &Y,
      Teuchos::ETransp mode,
      double alpha,
      double beta
      ) const;

  //virtual int
  //applyInverse(const Tpetra::MultiVector<double,int,int> &X,
  //             Tpetra::MultiVector<double,int,int> &Y
  //            ) const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>> getDomainMap() const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>> getRangeMap() const;

public:
  void
  rebuild(
      const std::map<std::string, double> & params,
      const Tpetra::Vector<double,int,int> &psi
      );

protected:
private:
  const std::shared_ptr<const Tpetra::Vector<double,int,int>>
  getAbsPsiSquared_(const Tpetra::Vector<double,int,int> &psi);

  void
  rebuildInverse_();

private:

  const std::shared_ptr<const nosh::mesh> mesh_;
  const std::shared_ptr<const nosh::scalar_field::base> thickness_;

  const std::shared_ptr<nosh::parameter_matrix::keo> regularizedkeo_;

  Teuchos::RCP<MueLu::TpetraOperator<double,int,int>> MueluPrec_;

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> timerRebuild0_;
  const Teuchos::RCP<Teuchos::Time> timerRebuild1_;
#endif

  const int num_cycles_;
};
} // namespace nosh

#endif // NOSH_KEOPRECONDITIONER_H
