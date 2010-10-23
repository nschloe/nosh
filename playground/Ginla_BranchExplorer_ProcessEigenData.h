/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef GINLA_BRANCHEXPLORER_PROCESSEIGENDATA_H
#define GINLA_BRANCHEXPLORER_PROCESSEIGENDATA_H

#include <LOCA_SaveEigenData_AbstractStrategy.H>
#include <Teuchos_Array.hpp>

// forward declarations
namespace Ginla {
  namespace IO {
    class StatsWriter;
    namespace EigenSaver {
      class Abstract;
    }
  }
}
namespace LOCA {
  class Stepper;
}
class Epetra_Vector;


namespace Ginla {
namespace BranchExplorer {

class ProcessEigenData:
        public LOCA::SaveEigenData::AbstractStrategy
{
public:

// Actually suggested interface:
//    EigenSaver(
//      const Teuchos::RCP<LOCA::GlobalData>& global_data,
//      const Teuchos::RCP<LOCA::P#ifndef GL_IO_SAVEEIGENDATA_Harameter::SublistParser>& topParams,
//      const Teuchos::RCP<Teuchos::ParameterList>& eigenParams      );

   // Constructor
   ProcessEigenData ( Teuchos::RCP<Teuchos::ParameterList>                      & eigenParamList,
                      const Teuchos::RCP<const Ginla::IO::EigenSaver::Abstract> & eigenSaver,
                      const Teuchos::RCP<Ginla::IO::StatsWriter>                & statsWriter
                    );
              
   virtual
   ~ProcessEigenData();

   virtual
   NOX::Abstract::Group::ReturnType
   save ( Teuchos::RCP<std::vector<double> >       & evals_r,
          Teuchos::RCP<std::vector<double> >       & evals_i,
          Teuchos::RCP<NOX::Abstract::MultiVector> & evecs_r,
          Teuchos::RCP<NOX::Abstract::MultiVector> & evecs_i
        );

   //! Returns a Teuchos::Array containing the eigenvectors of the eigenvalues
   //! that have just crossed the imaginary axis.
   Teuchos::Array<Epetra_Vector>
   getCrossingEigenvectors() const;
        
   void
   setLocaStepper( const Teuchos::RCP<LOCA::Stepper> locaStepper );

   // This function is necessary to break the circular dependency with the
   // LOCA_Stepper object to allow for a clean termination
   void
   releaseLocaStepper();

protected:
private:
    Teuchos::RCP<Teuchos::ParameterList> eigenParamList_;
    const Teuchos::RCP<const Ginla::IO::EigenSaver::Abstract> eigenSaver_;
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter_;
    Teuchos::RCP<LOCA::Stepper> locaStepper_;

    Teuchos::Array<Epetra_Vector> crossingEigenvectors_;

    unsigned int previousNumUnstableEigenvalues_;
    
    //! The minimum number of stable eigenvalues that is to be computed in each step.
    const unsigned int numComputeStableEigenvalues_;

    //! Maximum number of eigenvalues that are stored in \c eigenvaluesFilePath_.
    const unsigned int maxEigenvaluesSave_;
};

}

}

#endif // GINLA_BRANCHEXPLORER_PROCESSEIGENDATA_H
