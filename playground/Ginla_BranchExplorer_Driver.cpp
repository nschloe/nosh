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

#include "Ginla_BranchExplorer_Driver.h"

#include "Ginla_BranchExplorer_ProcessEigenData.h"

#include <Epetra_CrsMatrix.h>

#include <LOCA_Epetra_Factory.H>
#include <LOCA_Epetra_Group.H>
#include <LOCA_Stepper.H>
#include <LOCA_StatusTest_Combo.H>
#include <LOCA_StatusTest_MaxIters.H>

#include <NOX_StatusTest_MaxIters.H>
#include <NOX_StatusTest_NormF.H>
#include <NOX_StatusTest_Combo.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

#include "Ginla_Typedefs.h"
#include "Recti_Grid_Uniform.h"
#include "Ginla_MagneticVectorPotential_Centered.h"
#include "Ginla_Operator_BCCentral.h"
#include "Ginla_LocaSystem_Bordered.h"
#include "Ginla_Helpers.h"
#include "Ginla_StatusTest_Loop.h"
#include "Ginla_StatusTest_Turnaround.h"

#include "Ginla_IO_StatsWriter.h"
#include "Ginla_IO_StateWriter.h"

// ============================================================================
Ginla::BranchExplorer::Driver::
Driver():
   contFileBaseName_( "continuation" ),
   outputFormat_( "VTI" ),
   contDataFileName_( "continuationData.dat" )
{ 
}
// ============================================================================
Ginla::BranchExplorer::Driver::
~Driver()
{ 
}
// ============================================================================
void
Ginla::BranchExplorer::Driver::
run()
{
}
// ============================================================================
void
Ginla::BranchExplorer::Driver::
oneContinuation( const Teuchos::RCP<Epetra_Comm>            & eComm,
                 const std::string                          & outputDirectory,
                 const Teuchos::RCP<const ComplexVector>    & psi,
                 const Teuchos::RCP<Recti::Grid::Uniform>   & grid,
                 const Teuchos::ParameterList               & glParameters,
                 const Teuchos::RCP<Teuchos::ParameterList> & paramList
               )
{
    // ------------------------------------------------------------------------
    // create the problem
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( glParameters.get<double> ( "H0" ),
                                                                      glParameters.get<double> ( "scaling" ) ) );

    // create the operator
    Teuchos::RCP<Ginla::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::Operator::BCCentral ( grid, A, psi->getMap(), psi->getMap() ) );


    std::string fn = outputDirectory + "/" + contDataFileName_;
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter = 
        Teuchos::rcp( new Ginla::IO::StatsWriter( fn ) );

    Teuchos::RCP<Ginla::LocaSystem::Bordered> glsystem;

    const Teuchos::ParameterList & stepperList = paramList->sublist ( "LOCA" )
                                                           .sublist ( "Stepper" );
    int maxLocaSteps = stepperList.get<int> ( "Max Steps" );
    
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter = 
        Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory,
                                                  contFileBaseName_,
                                                  outputFormat_,
                                                  maxLocaSteps ) );

    glsystem = Teuchos::rcp ( new Ginla::LocaSystem::Bordered ( glOperator,
                                                                eComm,
                                                                psi->getMap(),
                                                                statsWriter,
                                                                stateWriter ) );
    // ------------------------------------------------------------------------
    // set the initial value from glParameters
    std::string contParam = stepperList.get<string> ( "Continuation Parameter" );
    TEST_FOR_EXCEPTION ( !glParameters.isParameter ( contParam ),
                         std::logic_error,
                         "Parameter \"" << contParam << "\" given as continuation parameter, but doesn't exist"
                         << "in the glParameters list." );

    // TODO Do we really not need this?
//     stepperList.set ( "Initial Value", glParameters.get<double> ( contParam ) );

    // ------------------------------------------------------------------------
    Teuchos::ParameterList outputList;
    outputList = paramList->sublist ( "Output", true );
    // ------------------------------------------------------------------------    
    // Create the necessary LOCA objects.
    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
        Teuchos::rcp ( new LOCA::Epetra::Factory );

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
        LOCA::createGlobalData ( paramList, epetraFactory );
    // ------------------------------------------------------------------------
    // Eigenstrategy
    Teuchos::ParameterList& eigenList =
        paramList->sublist ( "LOCA" ).sublist ( "Stepper" ) .sublist ( "Eigensolver" );
    Teuchos::RCP<Teuchos::ParameterList> eigenListPtr =
        Teuchos::rcpFromRef ( eigenList );
    std::string eigenvaluesFileName =
        outputDirectory  + "/" + outputList.get<string> ( "Eigenvalues file name" );
    std::string eigenstateFileNameAppendix =
        outputList.get<string> ( "Eigenstate file name appendix" );

    Teuchos::RCP<Ginla::IO::StatsWriter> eigenStatsWriter =
        Teuchos::rcp( new Ginla::IO::StatsWriter( eigenvaluesFileName ) );

    Teuchos::RCP<Ginla::BranchExplorer::ProcessEigenData> eigenProcessor =    
        Teuchos::rcp( new Ginla::BranchExplorer::ProcessEigenData ( eigenListPtr,
                                                                    glsystem,
                                                                    eigenStatsWriter ) );
                                                                                

    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy =
        eigenProcessor;
    eigenList.set ( "Save Eigen Data Method", "User-Defined" );
    eigenList.set ( "User-Defined Save Eigen Data Name", "glSaveEigenDataStrategy" );
    eigenList.set ( "glSaveEigenDataStrategy", glSaveEigenDataStrategy );
    // ------------------------------------------------------------------------
    // Create all possible Epetra_Operators.
    Teuchos::RCP<Epetra_RowMatrix> J = glsystem->getJacobian();
//  Teuchos::RCP<Epetra_RowMatrix> M = glsystem->getPreconditioner();

    // Create the linear system.
    // Use the TimeDependent interface for computation of shifted matrices.
    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = glsystem;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glsystem;
//  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = glsystem;

    Teuchos::ParameterList& nlPrintParams =
        paramList->sublist ( "NOX" ) .sublist ( "Printing" );

    Teuchos::ParameterList& lsParams =
        paramList->sublist ( "NOX" )
                  .sublist ( "Direction" )
                  .sublist ( "Newton" )
                  .sublist ( "Linear Solver" );

//  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(
//      new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams, iJac, J, iPrec, M, *soln));

    NOX::Epetra::Vector cloneVector( Epetra_Vector( *glsystem->getExtendedMap() ) );
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
        Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( nlPrintParams,
                                                              lsParams,
                                                              iReq,
                                                              iJac,
                                                              J,
                                                              cloneVector ) );

    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = glsystem;
    // ------------------------------------------------------------------------
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    // translate parameters into a LOCA list:
    LOCA::ParameterVector locaParams =
          *(Ginla::Helpers::teuchosParameterList2locaParameterVector( glParameters ));

    // get initial guess
    double chi = 0.0;
    NOX::Epetra::Vector initialGuess ( glsystem->createInterfaceState( psi, chi ),
                                       NOX::Epetra::Vector::CreateView
                                     );

    Teuchos::RCP<LOCA::Epetra::Group> grp =
        Teuchos::rcp ( new LOCA::Epetra::Group ( globalData,
                                                 nlPrintParams,
                                                 iTime,
                                                 initialGuess,
                                                 linSys,
                                                 linSys,
                                                 locaParams ) );

    grp->setParams ( locaParams );
    // ------------------------------------------------------------------------
    // Set up the NOX status tests
    Teuchos::ParameterList& noxList = paramList->sublist ( "NOX", true );
    double tol = noxList.get<double> ( "Tolerance" );
    int maxNonlinarSteps = noxList.get<int> ( "Max steps" );
    Teuchos::RCP<NOX::StatusTest::NormF> normF =
        Teuchos::rcp ( new NOX::StatusTest::NormF ( tol ) );
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters =
        Teuchos::rcp ( new NOX::StatusTest::MaxIters ( maxNonlinarSteps ) );
    Teuchos::RCP<NOX::StatusTest::Generic> comboOR =
        Teuchos::rcp ( new NOX::StatusTest::Combo ( NOX::StatusTest::Combo::OR,
                                                    normF,
                                                    maxIters ) );
    // ------------------------------------------------------------------------
    // Set up the LOCA status tests
    Teuchos::RCP<LOCA::StatusTest::Abstract> maxLocaStepsTest =
           Teuchos::rcp ( new LOCA::StatusTest::MaxIters ( maxLocaSteps ) );
    Teuchos::RCP<LOCA::StatusTest::Abstract> loopTest =
                     Teuchos::rcp ( new Ginla::StatusTest::Loop ( glsystem,
                                                                  grid ) );                                                                 
    Teuchos::RCP<LOCA::StatusTest::Abstract> turnaroundTest =
                     Teuchos::rcp ( new Ginla::StatusTest::Turnaround () );
    Teuchos::RCP<LOCA::StatusTest::Combo> locaCombo =
        Teuchos::rcp ( new LOCA::StatusTest::Combo ( LOCA::StatusTest::Combo::OR ) );
    locaCombo->addStatusTest( maxLocaStepsTest );
    locaCombo->addStatusTest( loopTest );
    locaCombo->addStatusTest( turnaroundTest );
    // ------------------------------------------------------------------------
    // Create the stepper
    Teuchos::RCP<LOCA::Stepper> stepper =
        Teuchos::rcp ( new LOCA::Stepper ( globalData,
                                           grp,
                                           locaCombo,
                                           comboOR,
                                           paramList ) );
    // ------------------------------------------------------------------------
    // pass pointer to stepper to glsystem to be able to read stats from the stepper in there
    glsystem->setLocaStepper ( stepper );
    eigenProcessor->setLocaStepper ( stepper );
    // ------------------------------------------------------------------------
    // Perform continuation run
    TEUCHOS_ASSERT_EQUALITY( 0, stepper->run() );
    // ------------------------------------------------------------------------
    // clean up
    LOCA::destroyGlobalData ( globalData );
    glsystem->releaseLocaStepper();
    eigenProcessor->releaseLocaStepper();
    // ------------------------------------------------------------------------
    
    return;
}
// ============================================================================