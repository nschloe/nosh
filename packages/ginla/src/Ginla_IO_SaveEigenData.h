#ifndef GINLA_IO_SAVEEIGENDATA_H
#define GINLA_IO_SAVEEIGENDATA_H

#include <Teuchos_RCP.hpp>

#include <LOCA_SaveEigenData_AbstractStrategy.H>
#include <LOCA_Parameter_SublistParser.H>
#include <LOCA_Stepper.H>
#include <NOX_Epetra_Vector.H>

#include <Teuchos_ParameterList.hpp>

#include "Komplex2_LinearProblem.h"

#include "Ginla_IO_EigenSaver_Abstract.h"
#include "Ginla_IO_StatsWriter.h"

// forward declarations
namespace Ginla {
  namespace IO {
    class StateWriter;
  }
  class StateTranslator;
}

namespace Ginla {
  namespace IO {

class SaveEigenData:
        public LOCA::SaveEigenData::AbstractStrategy
{

 public:

// Actually suggested interface:
//    EigenSaver(
//      const Teuchos::RCP<LOCA::GlobalData>& global_data,
//      const Teuchos::RCP<LOCA::P#ifndef GL_IO_SAVEEIGENDATA_Harameter::SublistParser>& topParams,
//      const Teuchos::RCP<Teuchos::ParameterList>& eigenParams      );

   // Constructor
   SaveEigenData ( Teuchos::ParameterList                           & eigenParamList,
                   const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator,
                   const Teuchos::RCP<const Ginla::IO::StateWriter> & stateWriter,
                   const Teuchos::RCP<Ginla::IO::StatsWriter>       & statsWriter
                 );
              
   virtual
   ~SaveEigenData();

   virtual
   NOX::Abstract::Group::ReturnType
   save ( Teuchos::RCP<std::vector<double> >       & evals_r,
          Teuchos::RCP<std::vector<double> >       & evals_i,
          Teuchos::RCP<NOX::Abstract::MultiVector> & evecs_r,
          Teuchos::RCP<NOX::Abstract::MultiVector> & evecs_i
        );

   void
   setLocaStepper( const Teuchos::RCP<LOCA::Stepper> locaStepper );

   // This function is necessary to break the circular dependency with the
   // LOCA_Stepper object to allow for a clean termination
   void
   releaseLocaStepper();

  protected:
  private:
    Teuchos::RCP<Teuchos::ParameterList> eigenParamListPtr_;
    const Teuchos::RCP<const Ginla::StateTranslator> stateTranslator_;
    const Teuchos::RCP<const Ginla::IO::StateWriter> stateWriter_;
    const Teuchos::RCP<const Komplex2::LinearProblem> komplex_;
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter_;
    Teuchos::RCP<LOCA::Stepper> locaStepper_;

    //! The minimum number of stable eigenvalues that is to be computed in each step.
    unsigned int numComputeStableEigenvalues_;

    //! Maximum number of eigenvalues that are stored in \c eigenvaluesFilePath_.
    unsigned int maxEigenvaluesSave_;
    unsigned int maxNumDigits_;
};

  } // namespace IO
} // namespace GL

#endif // GINLA_IO_SAVEEIGENDATA_H
