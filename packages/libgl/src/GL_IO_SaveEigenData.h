#ifndef GL_IO_SAVEEIGENDATA_H
#define GL_IO_SAVEEIGENDATA_H

#include <Teuchos_RCP.hpp>

#include <LOCA_SaveEigenData_AbstractStrategy.H>
#include <LOCA_Parameter_SublistParser.H>
#include <LOCA_Stepper.H>
#include <NOX_Epetra_Vector.H>

#include <Teuchos_ParameterList.hpp>

#include "AbstractStateWriter.h"

namespace GL {
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
   SaveEigenData( Teuchos::RCP<Teuchos::ParameterList> & eigenParams,
                  const std::string outputDir,
                  const std::string eigenvaluesFileName,
                  const std::string contFileBaseName,
                  const std::string eigenstateFileNameAppendix,
                  const Teuchos::RCP<AbstractStateWriter> glSys,
                  const unsigned int maxNumDigits
                );
               
   virtual
   ~SaveEigenData();

   virtual NOX::Abstract::Group::ReturnType
   save ( Teuchos::RCP<std::vector<double> >       &evals_r,
          Teuchos::RCP<std::vector<double> >       &evals_i,
          Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_r,
          Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_i  );

   void
   setLocaStepper( const Teuchos::RCP<LOCA::Stepper> locaStepper );

        // This function is necessary to break the circular dependency with the
        // LOCA_Stepper object to allow for a clean termination
        void
        releaseLocaStepper();

 protected:
  private:
    Teuchos::RCP<Teuchos::ParameterList> eigenParamList_;
    std::string outputDir_;
    std::string eigenvaluesFilePath_;
    std::string contFileBaseName_;
    std::string eigenstateFileNameAppendix_;
    Teuchos::RCP<AbstractStateWriter> stateWriter_;
    Teuchos::RCP<LOCA::Stepper> locaStepper_;

    //! The minimum number of stable eigenvalues that is to be computed in each step.
    unsigned int numComputeStableEigenvalues_;

    //! Maximum number of eigenvalues that are stored in \c eigenvaluesFilePath_.
    unsigned int maxEigenvaluesSave_;
    unsigned int maxNumDigits_;
};

  } // namespace IO
} // namespace GL

#endif // GL_IO_SAVEEIGENDATA_H
