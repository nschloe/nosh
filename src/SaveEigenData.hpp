// @HEADER
//
//    Helper class for computation and I/O of eigenvalues and -vectors.
//    Copyright (C) 2009--2012  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
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

#include "CsvWriter.hpp"
// =============================================================================
// forward declarations
namespace Nosh
{
namespace ModelEvaluator
{
class Virtual;
}
} // namespace Nosh
// =============================================================================
namespace Nosh
{

class SaveEigenData :
  public LOCA::SaveEigenData::AbstractStrategy
{
public:
// Actually suggested interface:
//    EigenSaver(
//      const std::shared_ptr<LOCA::GlobalData>& global_data,
//      const std::shared_ptr<LOCA::P#ifndef GL_SAVEEIGENDATA_Harameter::SublistParser>& topParams,
//      const std::shared_ptr<Teuchos::ParameterList>& eigenParams     );

  // Constructor
  SaveEigenData(
      Teuchos::ParameterList &eigenParamList,
      const std::shared_ptr<const Nosh::ModelEvaluator::Virtual> &modelEval,
      const std::string & fileName
      );

  virtual
  ~SaveEigenData();

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
  Teuchos::RCP<Teuchos::ParameterList> eigenParamListPtr_;
  const std::shared_ptr<const Nosh::ModelEvaluator::Virtual> modelEval_;
  Nosh::CsvWriter csvWriter_;
  std::shared_ptr<LOCA::Stepper> locaStepper_;

  //! The minimum number of stable eigenvalues that is to be computed in each
  //! step.
  unsigned int numComputeStableEigenvalues_;
};
} // namespace Nosh

#endif // NOSH_SAVEEIGENDATA_H
