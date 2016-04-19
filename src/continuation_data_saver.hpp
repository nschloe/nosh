#ifndef NOSH_CONTINUATION_DATA_SAVER
#define NOSH_CONTINUATION_DATA_SAVER

#include <LOCA_Thyra_SaveDataStrategy.H>
#include <NOX_Thyra_Vector.H>
#include <Thyra_TpetraThyraWrappers.hpp>

namespace nosh {
class continuation_data_saver: public LOCA::Thyra::SaveDataStrategy
{
  public:
  explicit continuation_data_saver(
      const std::shared_ptr<nosh::mesh> & mesh
      ):
    mesh_(mesh)
  {
  };

  virtual ~continuation_data_saver() {};

  virtual
  void
  saveSolution(
      const NOX::Abstract::Vector &x,
      double p
      )
  {
    // Get NOX::Thyra::Vector
    const auto x_nox_thyra = dynamic_cast<const NOX::Thyra::Vector*>(&x);
    TEUCHOS_ASSERT(x_nox_thyra != nullptr);
    const auto x_thyra = x_nox_thyra->getThyraRCPVector();
    auto x_tpetra =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getConstTpetraVector(
          x_thyra
          );
    std::cout << "XXX " << p << std::endl;
  }

  private:
  const std::shared_ptr<nosh::mesh> mesh_;
};
}  // namespace nosh
#endif  // NOSH_CONTINUATION_DATA_SAVER
