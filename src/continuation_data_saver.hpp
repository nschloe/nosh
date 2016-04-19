#ifndef NOSH_CONTINUATION_DATA_SAVER
#define NOSH_CONTINUATION_DATA_SAVER

#include <LOCA_Thyra_SaveDataStrategy.H>
#include <NOX_Thyra_Vector.H>
#include <Thyra_TpetraThyraWrappers.hpp>

#include "function.hpp"

namespace nosh {
class continuation_data_saver: public LOCA::Thyra::SaveDataStrategy
{
  public:
  explicit continuation_data_saver(
      const std::shared_ptr<nosh::mesh> & mesh
      ):
    mesh_(mesh),
    index_(0)
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
    (void) p;

    // extract Tpetra vector
    const auto x_nox_thyra = dynamic_cast<const NOX::Thyra::Vector*>(&x);
    TEUCHOS_ASSERT(x_nox_thyra != nullptr);
    const auto x_thyra = x_nox_thyra->getThyraRCPVector();
    auto x_tpetra =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getConstTpetraVector(
          x_thyra
          );

    std::ostringstream index_stream;
    index_stream << std::setw(4) << std::setfill('0') << index_;
    std::ostringstream filename;
    filename << "out" << index_stream.str() << ".h5m";

    nosh::write(Teuchos::get_shared_ptr(x_tpetra), mesh_, filename.str());

    index_++;
  }

  private:
  const std::shared_ptr<nosh::mesh> mesh_;
  size_t index_;
};
}  // namespace nosh
#endif  // NOSH_CONTINUATION_DATA_SAVER
