// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Epetra_LinearProblem.h>

#include <ml_epetra_preconditioner.h>

//#include "BelosConfigDefs.hpp"
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>

#include "Ginla_StkMeshReader.hpp"

#include "Ginla_State.hpp"
#include "Ginla_ModelEvaluator.hpp"
#include "Ginla_KeoContainer.hpp"
#include "Ginla_MagneticVectorPotential_ExplicitValues.hpp"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// =============================================================================
typedef double                           ST;
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
// =============================================================================
int main ( int argc, char *argv[] )
{
  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init( &argc, &argv );
#endif

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
      Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm>  eComm =
           Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

    const Teuchos::RCP<Teuchos::FancyOStream> out =
        Teuchos::VerboseObjectBase::getDefaultOStream();

    bool success = true;
    try
    {
      // ===========================================================================
      // handle command line arguments
      Teuchos::CommandLineProcessor My_CLP;

      My_CLP.setDocString (
              "Linear solver testbed for KEO and Jacobian operator.\n"
      );

      std::string inputFileName( "" );
      My_CLP.setOption ( "input", &inputFileName, "Input state file", true );

      bool verbose = true;
      My_CLP.setOption("verbose","quiet",&verbose,"Print messages and results.");

      bool isPrec = true;
      My_CLP.setOption("prec","noprec",&isPrec,"Use a preconditioner.");

      int frequency = 10;
      My_CLP.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");

      // print warning for unrecognized arguments
      My_CLP.recogniseAllOptions ( true );

      // finally, parse the command line
      TEUCHOS_ASSERT_EQUALITY( My_CLP.parse ( argc, argv ),
                               Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL
                             );
      // =========================================================================
      // Read the data from the file.
      Teuchos::ParameterList data;
      Teuchos::RCP<Teuchos::Time> readTime = Teuchos::TimeMonitor::getNewTimer("Data I/O");
      {
      Teuchos::TimeMonitor tm(*readTime);
      Ginla::StkMeshRead( *eComm, inputFileName, data );
      }

      // Cast the data into something more accessible.
      Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
      Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
      Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

      Teuchos::RCP<Ginla::MagneticVectorPotential::ExplicitValues> mvp;
      double mu;
      Teuchos::RCP<Teuchos::Time> mvpConstructTime = Teuchos::TimeMonitor::getNewTimer("MVP construction");
      {
          Teuchos::TimeMonitor tm(*mvpConstructTime);
          mu = problemParameters.get<double> ( "mu" );
          mu = 1.0e-2;
          mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential::ExplicitValues ( mesh, mvpValues, mu ) );
          //mvp->initializeEdgeMidpointProjectionCache_();
      }

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );

      Teuchos::RCP<Ginla::KeoContainer> keoContainer =
          Teuchos::rcp( new Ginla::KeoContainer( mesh, thickness, mvp ) );

      // create the kinetic energy operator
      Teuchos::RCP<Epetra_CrsMatrix> keoMatrix;
      Teuchos::RCP<Teuchos::Time> keoConstructTime = Teuchos::TimeMonitor::getNewTimer("Matrix construction");
      {
          Teuchos::TimeMonitor tm(*keoConstructTime);
          keoContainer->updateParameters( mvpParameters );
          // Copy over the matrix
          // TODO Well, we'll need to copy here as the KEO is otherwise negative
          // definite.
          keoMatrix = Teuchos::rcp( new Epetra_CrsMatrix( *keoContainer->getKeo() ) );
      }
      // Make sure the matrix is indeed positive definite, and not
      // negative definite. Belos needs that (2010-11-05).
      TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->Scale( -1.0 ) );

      // create initial guess and right-hand side
      bool zeroOut = true;
      Teuchos::RCP<Epetra_Vector> epetra_x =
              Teuchos::rcp( new Epetra_Vector( keoMatrix->OperatorDomainMap(), zeroOut ) );
      Teuchos::RCP<Epetra_MultiVector> epetra_b =
              Teuchos::rcp( new Epetra_Vector( keoMatrix->OperatorRangeMap(), 1 ) );
      //epetra_b->Random();
      TEUCHOS_ASSERT_EQUALITY( 0, epetra_b->PutScalar( 1.0 ) );

//      // check that the matrix is indeed symmetric
//      for ( int i=0; i<epetra_b->GlobalLength(); i++ )
//      {
//          int numEntries = 0;
//          double *values;
//          int *indices;
//          keoMatrix->ExtractMyRowView( i, numEntries, values, indices );
//          for ( int k=0; k<numEntries; k++ )
//          {
//              int j = indices[k];
//              if ( j>i )
//              {
//                  // get element [j,i]
//                  int numEntries2 = 0;
//                  double *values2;
//                  int *indices2;
//                  keoMatrix->ExtractMyRowView( j, numEntries2, values2, indices2 );
//                  // check if element [j,i] is there at all
//                  int kk = -1;
//                  for ( int k2=0; k2<numEntries2; k2++ )
//                  {
//                      if ( indices2[k2] == i )
//                      {
//                          kk = k2;
//                          break;
//                      }
//                  }
//                  // make sure it's been found
//                  TEUCHOS_ASSERT_INEQUALITY( kk, !=, -1 );
//                  // check that values are the same
//                  TEUCHOS_ASSERT_EQUALITY( values[k], values2[kk] );
//              }
//          }
//      }
//      *out << "Matrix appears to be symmetric" << std::endl;
      // -----------------------------------------------------------------------
      // Belos part
      Teuchos::ParameterList belosList;
      belosList.set( "Convergence Tolerance", 1.0e-15 );  // Relative convergence tolerance requested
      belosList.set( "Maximum Iterations", 10000 );
      if (verbose) {
        belosList.set( "Verbosity",
                       Belos::Errors +
                       Belos::Warnings +
                       Belos::IterationDetails +
                       Belos::FinalSummary +
                       Belos::Debug +
                       Belos::TimingDetails //+
//                       Belos::StatusTestDetails
                     );
        if (frequency > 0)
          belosList.set( "Output Frequency", frequency );
      }
      else
        belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

      // Only print on the zero processor
      const bool proc_verbose = verbose && (eComm->MyPID()==0);

      // Construct an unpreconditioned linear problem instance.
      Belos::LinearProblem<double,MV,OP> problem( keoMatrix, epetra_x, epetra_b );
      bool set = problem.setProblem();
      TEUCHOS_TEST_FOR_EXCEPTION( !set,
                                  std::runtime_error,
                                  "ERROR:  Belos::LinearProblem failed to set up correctly!" );
      // -----------------------------------------------------------------------
      // create preconditioner
      Teuchos::RCP<Teuchos::Time> precConstructTime = Teuchos::TimeMonitor::getNewTimer("Create preconditioner");
      if ( isPrec )
      {
          Teuchos::TimeMonitor tm(*precConstructTime);
          Teuchos::ParameterList MLList;
          ML_Epetra::SetDefaults( "SA", MLList );
          MLList.set("ML output", 0);
          MLList.set("max levels", 10);
          MLList.set("increasing or decreasing", "increasing");
          MLList.set("aggregation: type", "Uncoupled");
          MLList.set("smoother: type", "Chebyshev"); // "block Gauss-Seidel" "Chebyshev"
//           MLList.set("aggregation: threshold", 0.0);
          MLList.set("smoother: sweeps", 3);
          MLList.set("smoother: pre or post", "both");
          MLList.set("coarse: type", "Amesos-KLU");
          MLList.set("PDE equations", 2);
          Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec =
                      Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(*keoMatrix, MLList) );
          MLPrec->PrintUnused(0);

          Teuchos::RCP<Epetra_Operator> Prec =
                  Teuchos::rcp(  new ML_Epetra::MultiLevelPreconditioner(*keoMatrix, MLList) );
          TEUCHOS_ASSERT( !Prec.is_null() );

          // Create the Belos preconditioned operator from the preconditioner.
          // NOTE:  This is necessary because Belos expects an operator to apply the
          //        preconditioner with Apply() NOT ApplyInverse().
          Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( Prec ) );
          problem.setLeftPrec( belosPrec );
      }
      // -----------------------------------------------------------------------
      // Create an iterative solver manager.
      Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver
              = Teuchos::rcp( new Belos::PseudoBlockCGSolMgr<double,MV,OP>( Teuchos::rcp(&problem,false),
                                                                            Teuchos::rcp(&belosList,false)
                                                                          )
                   );
//       RCP< Belos::SolverManager<double,MV,OP> > newSolver
//               = rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>( rcp(&problem,false),
//                                                                       rcp(&belosList,false)
//                                                                     )
//                    );

      // Perform solve
      Belos::ReturnType ret = newSolver->solve();
      success = ret==Belos::Converged;
      // -----------------------------------------------------------------------
      //
      // Compute actual residuals.
      //
      bool badRes = false;
      std::vector<double> actual_resids( 1 );
      std::vector<double> rhs_norm( 1 );
      Epetra_Vector resid( keoMatrix->OperatorRangeMap() );
      OPT::Apply( *keoMatrix, *epetra_x, resid );
      MVT::MvAddMv( -1.0, resid, 1.0, *epetra_b, resid );
      MVT::MvNorm( resid, actual_resids );
      MVT::MvNorm( *epetra_b, rhs_norm );
      if (proc_verbose) {
        *out<< "---------- Actual Residuals (normalized) ----------" <<std::endl<<std::endl;
        for ( int i=0; i<1; i++) {
          double actRes = actual_resids[i]/rhs_norm[i];
          *out << "Problem " << i << " : \t" << actRes << std::endl;
          if (actRes > 1.0e-10) badRes = true;
        }
      }
      // -----------------------------------------------------------------------
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================

