/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Sch\"omer

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

#ifndef GINLA_LOCASYSTEM_VIRTUAL_H
#define GINLA_LOCASYSTEM_VIRTUAL_H


#include "Ginla_Komplex.h"
#include "Ginla_State.h"
#include "Ginla_IO_EigenSaver_Abstract.h"

#include <NOX_Epetra_Interface_Jacobian.H> // NOX base class
#include <NOX_Epetra_Interface_Preconditioner.H> // NOX base class
#include <LOCA_Epetra_Interface_TimeDependent.H> // LOCA base class

// // forward declarations
// namespace Ginla {
//   namespace Operator {
//     class Virtual;
//   }
//   namespace IO {
//     class StatsWriter;
//     class StateWriter;
//   }
//   namespace Perturbation {
//     class Virtual;
//   }
// }
namespace LOCA {
  class Stepper;
}



namespace Ginla {

namespace LocaSystem {

class Virtual:
            public Ginla::IO::EigenSaver::Abstract,
            public NOX::Epetra::Interface::Jacobian,
            public NOX::Epetra::Interface::Preconditioner,
            public LOCA::Epetra::Interface::TimeDependent
{
public:

    //! Constructor with initial guess.
    Virtual ();

    //! Destructor
    ~Virtual();

    //! Evaluate the Ginzburg--Landau functions at a given state defined
    //! by the input vector x.
    virtual bool
    computeF ( const Epetra_Vector & x,
               Epetra_Vector       & F,
               const NOX::Epetra::Interface::Required::FillType fillFlag = Residual ) = 0;
               
    //! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
    //! at a given state defined by the input vector x.
    virtual bool
    computeJacobian ( const Epetra_Vector & x,
                      Epetra_Operator     & Jac ) = 0;

    virtual bool
    computeShiftedMatrix ( double alpha,
                           double beta,
                           const Epetra_Vector & x,
                           Epetra_Operator     & A ) = 0;

    //! Dummy preconditioner function. So far does nothing but throwing
    //! an exception when called.
    virtual bool
    computePreconditioner ( const Epetra_Vector & x,
                            Epetra_Operator     & Prec,
                            Teuchos::ParameterList* precParams = 0 ) = 0;

    //! Returns the current Jacobian.
    //! @return Reference-counted pointer to the Jacobian.
    virtual Teuchos::RCP<Epetra_CrsMatrix>
    getJacobian() const = 0;

    virtual Teuchos::RCP<Epetra_CrsMatrix>
    getPreconditioner() const = 0;

    //! Set the problem parameters.
    virtual void
    setParameters ( const LOCA::ParameterVector &p ) = 0;

    //! Set directory to where all output gets written.
    //! @param directory Name of the directory.
    virtual void
    setOutputDir ( const string &directory ) = 0;

    //! Print the solution x along with the continuation parameter conParam
    //! to a file. This function is called internally by LOCA to print
    //! solutions of each continuation step.
    virtual void
    printSolution ( const Epetra_Vector & x,
                    double conParam ) = 0;
                    
    //! Used to print eigenvectors.
    virtual void
    printSolution ( const Epetra_Vector & x,
                    const std::string   & filenameAppendix
                  ) const = 0;

    virtual void
    setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper ) = 0;

    // This function is necessary to break the circular dependency with the
    // LOCA_Stepper object to allow for a clean termination
    virtual void
    releaseLocaStepper() = 0;

    virtual Teuchos::RCP<const Ginla::Komplex>
    getKomplex() const = 0;
    
    virtual Teuchos::RCP<const Epetra_Map>
    getMap() const = 0;
    
    virtual Teuchos::RCP<Epetra_Vector>
    createSystemVector( const Teuchos::ParameterList & p ) = 0;
    
//     Teuchos::RCP<ComplexVector>
//     extractPsi( const Epetra_Vector & x ) const;

    virtual Teuchos::RCP<Ginla::State>
    createState(  const Epetra_Vector & x ) const = 0;;

  protected:
  private:
};

}

}

#endif // GINLA_LOCASYSTEM_VIRTUAL_H