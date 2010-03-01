/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2009--2010 Nico Schl\"omer

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

#ifndef GLSYSTEMLOCACONSTRAINT_H
#define GLSYSTEMLOCACONSTRAINT_H

#include <LOCA_MultiContinuation_ConstraintInterfaceMVDX.H>

#include "glSystem.h"

class GlSystemLocaConstraint:
            public AbstractStateWriter,
            public NOX::Epetra::Interface::Jacobian,
            public NOX::Epetra::Interface::Preconditioner,
            public LOCA::Epetra::Interface::TimeDependent
//             public LOCA::MultiContinuation::ConstraintInterfaceMVDX
{
  public:

    //! Constructor with initial guess.
    GlSystemLocaConstraint ( GinzburgLandau::GinzburgLandau &gl,
                             const Teuchos::RCP<const Epetra_Comm> eComm,
                             const Teuchos::RCP<const ComplexVector> psi,
                             const std::string outputDir = "data",
                             const std::string outputDataFileName = "continuationData.dat",
                             const std::string outputFileFormat = "VTK",
                             const std::string solutionFileNameBase = "solutionStep",
                             const std::string nullvectorFileNameBase = "nullvectorStep",
                             const unsigned int maxStepNumberDecimals = 4 );
  
    //! Destructor.
    ~GlSystemLocaConstraint();

    //! Evaluate the Ginzburg--Landau functions at a given state defined
    //! by the input vector x.
    virtual bool
    computeF ( const Epetra_Vector &x,
               Epetra_Vector &F,
               const NOX::Epetra::Interface::Required::FillType fillFlag = Residual );

    //! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
    //! at a given state defined by the input vector x.
    virtual bool
    computeJacobian ( const Epetra_Vector &x, Epetra_Operator &Jac );

    virtual bool
    computeShiftedMatrix ( double alpha,
                           double beta,
                           const Epetra_Vector &x,
                           Epetra_Operator &A );

    //! Dummy preconditioner function. So far does nothing but throwing
    //! an exception when called.
    virtual bool
    computePreconditioner ( const Epetra_Vector    & x,
                            Epetra_Operator        & Prec,
                            Teuchos::ParameterList * precParams = 0 );

    //! Returns the current state. Not necessarily a solution to the problem.
    //! @return Reference-counted pointer to the current state.
    Teuchos::RCP<Epetra_Vector>
    getSolution() const;

    //! Returns the current Jacobian.
    //! @return Reference-counted pointer to the Jacobian.
    Teuchos::RCP<Epetra_CrsMatrix>
    getJacobian() const;

    Teuchos::RCP<Epetra_CrsMatrix>
    getPreconditioner() const;

    //! Set the problem parameters.
    virtual void
    setParameters ( const LOCA::ParameterVector & p );

    //! Set directory to where all output gets written.
    //! @param directory Name of the directory.
    void
    setOutputDir ( const string & directory );

    //! Print the solution x along with the continuation parameter conParam
    //! to a file. This function is called internally by LOCA to print
    //! solutions of each continuation step.
    virtual void
    printSolution ( const Epetra_Vector &x, double conParam );

    void
    setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper );

    // This function is necessary to break the circular dependency with the
    // LOCA_Stepper object to allow for a clean termination
    void
    releaseLocaStepper();

    //! Explicitly print the solution x along with the problem parameters
    //! to the file fileName.
    void
    writeSolutionToFile ( const Epetra_Vector & x,
                          const std::string   & filePath ) const;

    void
    writeAbstractStateToFile ( const Epetra_Vector & x,
                               const std::string   & filePath ) const;

    // TODO delete?
    const Teuchos::RCP<const GlKomplex>
    getGlKomplex() const;

    // TODO delete
    double
    getH0() const;
    void
    setH0 ( const double h0 );
    void
    setScaling ( const double h0 );
    void
    setChi ( const double h0 );

private:

    //! Print method for the continuation in one parameter.
    void
    printSolutionOneParameterContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                          ) const;

    //! Print method for turning point continuation continuation.
    void
    printSolutionTurningPointContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                          ) const;


    //! Write statistics about the current continuation step to the file
    //! \c outputDataFileName_ .
    void
    writeContinuationStats ( const int conStep,
                             const Teuchos::RCP<const ComplexVector> psi ) const;

private:

    GlSystem::GlSystem glSystem_;
};

#endif // GLSYSTEMLOCACONSTRAINT_H



