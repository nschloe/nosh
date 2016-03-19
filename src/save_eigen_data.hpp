#ifndef NOSH_SAVEEIGENDATA_H
#define NOSH_SAVEEIGENDATA_H
// =============================================================================
#include <string>
#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <LOCA_SaveEigenData_AbstractStrategy.H>
#include <LOCA_Parameter_SublistParser.H>
#include <LOCA_Stepper.H>

#include "csv_writer.hpp"
// =============================================================================
// forward declarations
namespace nosh
{
namespace model_evaluator
{
class base;
}
} // namespace nosh
// =============================================================================
namespace nosh
{

class save_eigen_data :
  public LOCA::SaveEigenData::AbstractStrategy
{
public:
// Actually suggested interface:
//    EigenSaver(
//      const std::shared_ptr<LOCA::Global_data>& global_data,
//      const std::shared_ptr<LOCA::P#ifndef GL_SAVEEIGENDATA_Harameter::SublistParser>& topParams,
//      const std::shared_ptr<Teuchos::ParameterList>& eigenParams     );

  // Constructor
  save_eigen_data(
      Teuchos::ParameterList &eigen_param_list,
      const std::shared_ptr<const nosh::model_evaluator::base> &model_eval,
      const std::string & file_name
      );

  virtual
  ~save_eigen_data();

  virtual
  NOX::Abstract::Group::ReturnType
  save(
      Teuchos::RCP<std::vector<double> > &evals_r,
      Teuchos::RCP<std::vector<double> > &evals_i,
      Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_r,
      Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_i
      );

  void
  setLocaStepper(const std::shared_ptr<LOCA::Stepper> locaStepper);

  // This function is necessary to break the circular dependency with the
  // LOCA_Stepper object to allow for a clean termination
  void
  releaseLocaStepper();

protected:
private:
  Teuchos::RCP<Teuchos::ParameterList> eigen_param_listPtr_;
  const std::shared_ptr<const nosh::model_evaluator::base> model_eval_;
  nosh::csv_writer csv_writer_;
  std::shared_ptr<LOCA::Stepper> locaStepper_;
};
} // namespace nosh

#endif // NOSH_SAVEEIGENDATA_H
