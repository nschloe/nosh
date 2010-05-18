/*
 * glNoxHelpers.cpp
 *
 *  Created on: Jan 17, 2010
 *      Author: Nico Schloemer
 */

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Comm.h>

#include <Tpetra_Vector.hpp>

#include <NOX.H>
#include <NOX_Epetra.H>

#include <NOX.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_Epetra.H>

// for the eigenvalue computation:
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBasicOutputManager.hpp>
#include <AnasaziBlockDavidsonSolMgr.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziLOBPCGSolMgr.hpp>
#include <AnasaziRTRSolMgr.hpp>
#include <AnasaziConfigDefs.hpp>
#include <AnasaziEpetraAdapter.hpp>
#include <AnasaziBasicSort.hpp>

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Ginla_IO_SaveNewtonData.h"
#include "Ginla_IO_StateWriter.h"

#include "Recti_Grid_Reader.h"
#include "Recti_Grid_UniformSquare.h"
#include "Recti_Grid_Square.h"

#include "Ginla_LocaSystem_Virtual.h"
#include "Ginla_LocaSystem_Default.h"
#include "Ginla_LocaSystem_Bordered.h"

#include "Ginla_MagneticVectorPotential_Centered.h"

#include "Recti_Domain_Factory.h"

#include "Ginla_IO_StatsWriter.h"

#include "Ginla_Operator_BCInner.h"
#include "Ginla_Operator_BCOuter.h"
#include "Ginla_Operator_BCCentral.h"

typedef Tpetra::Vector<std::complex<double>, Thyra::Ordinal> ComplexVector;

namespace glNoxHelpers
{
// =========================================================================
void
createGlSystem ( const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                 const Teuchos::RCP<const Epetra_Comm>         & eComm,
                 const std::string                             & fileName,
                 Teuchos::ParameterList                        & problemParameters,
                 Teuchos::RCP<Ginla::LocaSystem::Virtual>      & glSystem,
                 Teuchos::RCP<ComplexVector>                   & initialPsi,
                 Teuchos::RCP<Recti::Grid::Uniform>            & grid
               )
{
    Teuchos::ParameterList glParameters;

    Teuchos::RCP<ComplexMultiVector> psi;

    Recti::Grid::Reader::read ( comm, fileName, psi, grid, problemParameters );

    double h0      = problemParameters.get<double> ( "H0" );
    double scaling = problemParameters.get<double> ( "scaling" );

    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( h0, scaling ) );

    // create the operator
    Teuchos::RCP<Ginla::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::Operator::BCCentral ( grid, A, psi->getMap(), psi->getMap() ) );

    TEUCHOS_ASSERT_EQUALITY ( psi->getNumVectors(), 1 );
    initialPsi = psi->getVectorNonConst( 0 );
    
    std::string outputDirectory = "";
    std::string contDataFileName = "continuationData.dat";
    std::string contFileBaseName = "newtonStep";
    std::string outputFormat = "VTI";
    unsigned int maxIndex = 1000; // TODO replace by maxnumsteps
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter = Teuchos::null;
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter = 
        Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory,
                                                  contFileBaseName,
                                                  outputFormat,
                                                  maxIndex ) );
    
    glSystem = Teuchos::rcp ( new Ginla::LocaSystem::Default ( glOperator,
                                                               eComm,
                                                               psi->getMap(),
                                                               statsWriter,
                                                               stateWriter ) );

    return;
}
// =============================================================================
void
createGlSystem ( const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                 const Teuchos::RCP<const Epetra_Comm>         & eComm,
                 const unsigned int Nx,
                 const double scaling,
                 const double H0,
                 const Teuchos::ParameterList                & domainParameters,
                 Teuchos::ParameterList                      & problemParameters,
                 Teuchos::RCP<Ginla::LocaSystem::Virtual>    & glSystem,
                 Teuchos::RCP<ComplexVector>                 & initialPsi,
                 Teuchos::RCP<Recti::Grid::Uniform>          & grid
               )
{ 
    problemParameters.set ( "scaling", scaling );
    problemParameters.set ( "H0"     , H0 );
    problemParameters.set ( "chi"    , 0.0 );

    // create the domain
    Teuchos::RCP<Recti::Domain::Abstract> domain =
            Recti::Domain::buildDomain( domainParameters );
            
    // TODO Create GridConstructor with Nx
    // create the grid
    double h = 1.0 / Nx;
    grid = Teuchos::rcp ( new Recti::Grid::Uniform ( domain, h ) );
    
    grid->updateScaling ( scaling );
    
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( H0, scaling ) );

    // create an initial guess
    int numGlobalElements = grid->getNumGridPoints();
    int indexBase = 0;
    Teuchos::RCP<Tpetra::Map<Thyra::Ordinal> > map =
        Teuchos::rcp( new Tpetra::Map<Thyra::Ordinal>(numGlobalElements, indexBase, comm ) );

    // create the operator
    Teuchos::RCP<Ginla::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::Operator::BCCentral ( grid, A, map, map ) );
        
    initialPsi = Teuchos::rcp( new ComplexVector(map) );
    initialPsi->putScalar( double_complex(0.5,0.0) );
    
    std::string outputDirectory = "";
    std::string contDataFileName = "continuationData.dat";
    std::string contFileBaseName = "continuationStep";
    std::string outputFormat = "VTI";
    unsigned int maxIndex = 1000;

    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter = Teuchos::null;
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter = 
        Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory,
                                                  contFileBaseName,
                                                  outputFormat,
                                                  maxIndex ) );
            
    glSystem = Teuchos::rcp ( new Ginla::LocaSystem::Default ( glOperator,
                                                                eComm,
                                                                initialPsi->getMap(),
                                                                statsWriter,
                                                                stateWriter ) );
    return;
}
// =========================================================================
Teuchos::RCP<NOX::Epetra::Group>
createSolverGroup ( const Teuchos::RCP<Ginla::LocaSystem::Virtual>  & glSystem,
                    const Teuchos::RCP<Teuchos::ParameterList>      & nlParamsPtr,
                    const Teuchos::RCP<ComplexVector>               & initialPsi
                  )
{
    // Create all possible Epetra_Operators.
    Teuchos::RCP<Epetra_RowMatrix> Analytic = glSystem->getJacobian();

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = glSystem;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glSystem;

    Teuchos::ParameterList & printParams = nlParamsPtr->sublist ( "Printing" );

    Teuchos::ParameterList & lsParams    = nlParamsPtr->sublist ( "Direction" )
                                           .sublist ( "Newton" )
                                           .sublist ( "Linear Solver" );

    // Get initial solution
    Epetra_Vector cloneVector( *glSystem->getMap() );
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
        Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( printParams,
                                                              lsParams,
                                                              iReq,
                                                              iJac,
                                                              Analytic,
                                                              cloneVector ) );


    // get initial guess
    Teuchos::ParameterList p;
    p.set( "psi", initialPsi );
    p.set( "chi", 0.0 );
    NOX::Epetra::Vector initialGuess ( glSystem->createSystemVector( p ),
                                       NOX::Epetra::Vector::CreateView
                                     );

    // Create the Group
    return Teuchos::rcp ( new NOX::Epetra::Group ( printParams,
                                                   iReq,
                                                   initialGuess,
                                                   linSys ) );
}
// =========================================================================
Teuchos::RCP<NOX::StatusTest::Generic>
createConvergenceTest ( Teuchos::ParameterList & noxStatusList,
                        Teuchos::ParameterList & nlParamsPtr )
{
    NOX::StatusTest::Factory statusTestFactory;

    Teuchos::ParameterList & printParams = nlParamsPtr.sublist ( "Printing" );
    NOX::Utils outputUtils ( printParams );

    return statusTestFactory.buildStatusTests ( noxStatusList, outputUtils ) ;
}
// =============================================================================
// Create the solver
Teuchos::RCP<NOX::Solver::Generic>
createSolver ( const Teuchos::RCP<NOX::Epetra::Group>       grpPtr,
               const Teuchos::RCP<NOX::StatusTest::Generic> statusTest,
               const Teuchos::RCP<Teuchos::ParameterList>   nlParamsPtr )
{
    TEST_FOR_EXCEPTION ( grpPtr.is_null(),
                         std::logic_error,
                         "Group not initialized." );

    TEST_FOR_EXCEPTION ( statusTest.is_null(),
                         std::logic_error,
                         "Status test not initialized." );

    TEST_FOR_EXCEPTION ( nlParamsPtr.is_null(),
                         std::logic_error,
                         "Nonlinear solver parameters not initialized." );

    return NOX::Solver::buildSolver ( grpPtr,
                                      statusTest,
                                      nlParamsPtr );
}
// =========================================================================
double
computeJacobianConditionNumber ( const Teuchos::RCP<const NOX::Solver::Generic> solver,
                                 const Teuchos::RCP<      NOX::Epetra::Group>   grpPtr )
{
    double kappa;
    try
    {
        const NOX::Epetra::Group& finalGroup =
            dynamic_cast<const NOX::Epetra::Group&> ( solver->getSolutionGroup() );
        grpPtr->computeJacobian();
        grpPtr->computeJacobianConditionNumber ( 2000, 1e-2, 30, true );
        kappa = finalGroup.getJacobianConditionNumber();
    }
    catch ( std::exception& e )
    {
        std::cerr << e.what() << std::endl;
    }

    return kappa;
}
// =============================================================================
// compute some eigenvalues using anasazi
void
computeJacobianEigenvalues ( const Teuchos::RCP<const NOX::Solver::Generic> solver,
                             const Teuchos::RCP<      NOX::Epetra::Group>   grpPtr,
                             const int MyPID )
{
    bool debug = true; // even more increased verbosity

    std::string which ( "LR" ); // compute the rightmost eigenvalues

    bool boolret;

    typedef double ScalarType;
    typedef Teuchos::ScalarTraits<ScalarType>          SCT;
    typedef SCT::magnitudeType               MagnitudeType;
    typedef Epetra_MultiVector                          MV;
    typedef Epetra_Operator                             OP;
    typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
    typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

    const NOX::Epetra::Group& finalGroup =
        dynamic_cast<const NOX::Epetra::Group&> ( solver->getSolutionGroup() );

    const Epetra_Vector& finalSolution =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( finalGroup.getX() ) ).
        getEpetraVector();

    Epetra_BlockMap Map = finalSolution.Map();

    // shortname for the final Jacobian
    Teuchos::RCP<Epetra_Operator> J =
        grpPtr->getLinearSystem()->getJacobianOperator();

    // - - - - - - - - - - - - - - - - -
    // Start the block Arnoldi iteration
    // - - - - - - - - - - - - - - - - -
    // Variables used for the Block Krylov Schur Method
    int nev = 10;
    int blockSize = 1;
    int numBlocks = 30;
    int maxRestarts = 250;
    double tol = 1e-10;

    // Create a sort manager to pass into the block Krylov-Schur solver manager
    // -->  Make sure the reference-counted pointer is of type Anasazi::SortManager<>
    // -->  The block Krylov-Schur solver manager uses Anasazi::BasicSort<> by default,
    //      so you can also pass in the parameter "Which", instead of a sort manager.
    Teuchos::RCP<Anasazi::SortManager<MagnitudeType> > MySort =
        Teuchos::rcp ( new Anasazi::BasicSort<MagnitudeType> ( which ) );

    // Set verbosity level
    int verbosity = Anasazi::Errors + Anasazi::Warnings;
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
    if ( debug )
    {
        verbosity += Anasazi::Debug;
    }

    // Create parameter list to pass into solver manager
    Teuchos::ParameterList MyPL;
    MyPL.set ( "Verbosity", verbosity );
    MyPL.set ( "Sort Manager", MySort );
    //MyPL.set( "Which", which );
    MyPL.set ( "Block Size", blockSize );
    MyPL.set ( "Num Blocks", numBlocks );
    MyPL.set ( "Maximum Restarts", maxRestarts );
    //MyPL.set( "Step Size", stepSize );
    MyPL.set ( "Convergence Tolerance", tol );

    // Create an Epetra_MultiVector for an initial vector to start the solver.
    // Note:  This needs to have the same number of columns as the blocksize.
    Teuchos::RCP<Epetra_MultiVector> ivec =
        Teuchos::rcp ( new Epetra_MultiVector ( Map, blockSize ) );
    ivec->Random();

    // Create the eigenproblem.
    Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
        Teuchos::rcp ( new Anasazi::BasicEigenproblem<double, MV, OP> ( J, ivec ) );

    // Set the number of eigenvalues requested
    MyProblem->setNEV ( nev );

    // Inform the eigenproblem that you are finishing passing it information
    boolret = MyProblem->setProblem();
    TEST_FOR_EXCEPTION ( !boolret,
                         std::runtime_error,
                         "Anasazi::BasicEigenproblem::setProblem() returned with error." );

    // Initialize the Block Arnoldi solver
//       Anasazi::BlockDavidsonSolMgr<double, MV, OP>
//       Anasazi::LOBPCGSolMgr<double, MV, OP>
//       Anasazi::RTRSolMgr<double, MV, OP>
    Anasazi::BlockKrylovSchurSolMgr<double, MV, OP>
    MySolverMgr ( MyProblem, MyPL );

    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMgr.solve();
    if ( returnCode != Anasazi::Converged && MyPID==0 )
        cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;

    // Get the Ritz values from the eigensolver
//       std::vector<Anasazi::Value<double> > ritzValues = MySolverMgr.getRitzValues();

//       // Output computed eigenvalues and their direct residuals
//       if (MyPID==0) {
//         int numritz = (int)ritzValues.size();
//         cout.setf(std::ios_base::right, std::ios_base::adjustfield);
//         cout<<endl<< "Computed Ritz Values"<< endl;
//         if (MyProblem->isHermitian()) {
//           cout<< std::setw(16) << "Real Part"
//             << endl;
//           cout<<"-----------------------------------------------------------"<<endl;
//           for (int i=0; i<numritz; i++) {
//             cout<< std::setw(16) << ritzValues[i].realpart
//               << endl;
//           }
//           cout<<"-----------------------------------------------------------"<<endl;
//         }
//         else {
//           cout<< std::setw(16) << "Real Part"
//             << std::setw(16) << "Imag Part"
//             << endl;
//           cout<<"-----------------------------------------------------------"<<endl;
//           for (int i=0; i<numritz; i++) {
//             cout<< std::setw(16) << ritzValues[i].realpart
//               << std::setw(16) << ritzValues[i].imagpart
//               << endl;
//           }
//           cout<<"-----------------------------------------------------------"<<endl;
//           }
//         }

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
    std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
    Teuchos::RCP<MV> evecs = sol.Evecs;
    std::vector<int> index = sol.index;
//     int numev = sol.numVecs;

    Teuchos::RCP<Epetra_Vector> ev1 ( ( *evecs ) ( 0 ) );

//   glsystem->solutionToFile( *ev1,
//                             problemParameters,
//                             "data/ev1.vtk" );
//
//       cout<< std::setw(16) << "Real Part"
//           << std::setw(16) << "Imag Part" << endl;
//       cout<<"-----------------------------------------------------------"<<endl;
//       for (int i=0; i<numev; i++) {
//           cout<< std::setw(16) << evals[i].realpart
//             << std::setw(16) << evals[i].imagpart << endl;
//       }
    // -----------------------------------------------------------------------
}
// =========================================================================
void
printSolutionToFile ( const std::string                                     & outputDir,
                      const std::string                                     & fileBaseName,
                      const std::string                                     & outputFormat,
                      const Teuchos::RCP<const NOX::Solver::Generic>        & solver,
                      const Teuchos::RCP<const Ginla::LocaSystem::Virtual>  & glSystem
                    )                    
{
    const NOX::Epetra::Group & finalGroup =
        dynamic_cast<const NOX::Epetra::Group&> ( solver->getSolutionGroup() );

    const Epetra_Vector & finalSolution =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( finalGroup.getX() ) ).getEpetraVector();

    Teuchos::RCP<Ginla::State> state = glSystem->createState( finalSolution );
        
    unsigned int maxIndex = 0;
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter
        = Teuchos::rcp( new Ginla::IO::StateWriter( outputDir,
                                                    fileBaseName,
                                                    outputFormat ,
                                                    maxIndex ) );

    stateWriter->write( state, 0 );

    return;
}
// =========================================================================
int
checkConvergence ( const Teuchos::RCP<const NOX::Solver::Generic> solver )
{
    int status = 0;

//   // 1. Convergence
//   if ( solvStatus != NOX::StatusTest::Converged ) {
//       status = 1;
//       if ( printing.isPrintType ( NOX::Utils::Error ) )
//         printing.out() << "Nonlinear solver failed to converge!" << endl;
//   }

#ifndef HAVE_MPI
    // 2. Linear solve iterations (53) - SERIAL TEST ONLY!
    //    The number of linear iterations changes with # of procs.
    if ( const_cast<Teuchos::ParameterList&> ( solver->getList() )
         .sublist ( "Direction" )
         .sublist ( "Newton" )
         .sublist ( "Linear Solver" )
         .sublist ( "Output" )
         .get ( "Total Number of Linear Iterations",0 )
         != 53 )
    {
        status = 2;
    }
#endif

//  // 3. Nonlinear solve iterations (10)
//  if ( const_cast<Teuchos::ParameterList&> ( solver_->getList() )
//                                             .sublist ( "Output" )
//                                             .get ( "Nonlinear Iterations", 0 )
//        == maxNonlinearIterations_ )
//    status = 3;

    return status;
}
// =========================================================================
} // namespace glNoxHelpers
