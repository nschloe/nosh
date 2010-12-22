#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include <boost/filesystem.hpp>

#include "Ginla_EpetraFVM_StkMeshReader.h"

#include "Ginla_EpetraFVM_State.h"
#include "Ginla_EpetraFVM_ModelEvaluator.h"
#include "Ginla_EpetraFVM_KeoPreconditioner.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"
#include "Ginla_IO_NoxObserver.h"
#include "Ginla_IO_SaveEigenData.h"
#include "Ginla_MagneticVectorPotential_Spherical.h"
#include "Ginla_MagneticVectorPotential_Z.h"
#include "Ginla_MagneticVectorPotential_MagneticDot.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

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

    int status = 0;

    try
    {
      // ===========================================================================
      // handle command line arguments
      Teuchos::CommandLineProcessor My_CLP;

      My_CLP.setDocString (
          "Linear solver testbed for KEO and Jacobian operator.\n"
      );

      std::string inputFileName = "";
      My_CLP.setOption ( "input", &inputFileName,
                         "Input state file", true
                       );

      // print warning for unrecognized arguments
      My_CLP.recogniseAllOptions ( true );

      // do throw exceptions
      My_CLP.throwExceptions ( true );

      // finally, parse the command line
      My_CLP.parse ( argc, argv );
      // =========================================================================

      Teuchos::ParameterList              problemParameters;
      Teuchos::RCP<Epetra_Vector>         z = Teuchos::null;
      Teuchos::RCP<Ginla::EpetraFVM::StkMesh> mesh = Teuchos::null;

      Ginla::EpetraFVM::StkMeshRead( *eComm,
                                      inputFileName,
                                      z,
                                      mesh,
                                      problemParameters
                                    );
      double mu = problemParameters.get<double> ( "mu" );

      Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp =
              Teuchos::rcp ( new Ginla::MagneticVectorPotential::Z ( mesh, mu ) );

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );
      Teuchos::Tuple<double,3> scaling = Teuchos::tuple( 1.0, 1.0, 1.0 );
      double temperature = 1.0;

      // create the preconditioner
      Teuchos::RCP<Ginla::EpetraFVM::KeoPreconditioner> keoPrec =
              Teuchos::rcp( new  Ginla::EpetraFVM::KeoPreconditioner( mesh, mvp ) );
      keoPrec->rebuild( mvpParameters, scaling );

      // create the jacobian
      Teuchos::RCP<Ginla::EpetraFVM::JacobianOperator> jacobian =
              Teuchos::rcp( new Ginla::EpetraFVM::JacobianOperator( mesh, mvp ) );
      jacobian->rebuild( mvpParameters, scaling, temperature, z );

      // create initial guess and right-hand side
      Teuchos::RCP<Epetra_Vector> epetra_x = Teuchos::rcp( new Epetra_Vector( jacobian->OperatorDomainMap() ) );
      Teuchos::RCP<Epetra_Vector> epetra_b = Teuchos::rcp( new Epetra_Vector( jacobian->OperatorRangeMap() ) );
      epetra_b->PutScalar( 1.0 );

      // -----------------------------------------------------------------------
      // build the AztecOO problem
      Epetra_LinearProblem problem( &*jacobian, &*epetra_x, &*epetra_b );
      // make sure the problem is symmetric
      problem.AssertSymmetric();

      AztecOO solver( problem );
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      solver.SetAztecOption(AZ_precond, AZ_none);
//       solver.SetPrecOperator( &*keoPrec );

      //solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
      solver.SetAztecOption(AZ_solver, AZ_gmres);
      //solver.SetAztecOption(AZ_solver, AZ_cg);
//       solver.SetAztecOption(AZ_solver, AZ_bicgstab);
      //solver.SetAztecOption(AZ_scaling, 8);
      //solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
      //solver.SetAztecOption(AZ_output, 1);
      //solver.SetAztecOption(AZ_reorder, 0);
      //solver.SetAztecOption(AZ_graph_fill, 3);
      //solver.SetAztecOption(AZ_overlap, 0);
      //solver.SetAztecOption(AZ_poly_ord, 9);
      //solver.SetAztecParam(AZ_ilut_fill, 4.0);
      //solver.SetAztecParam(AZ_drop, 0.0);
      solver.SetAztecOption(AZ_output, 1);
      //double rthresh = 1.4;
      //cout << "Rel threshold = " << rthresh << endl;
      //solver.SetAztecParam(AZ_rthresh, rthresh);
      //double athresh = 10.0;
      //cout << "Abs threshold = " << athresh << endl;
      //solver.SetAztecParam(AZ_athresh, athresh);
      //solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);

      int Niters = 10000;
      //solver.SetAztecOption(AZ_kspace, Niters);

      // do the iteration
      solver.Iterate(Niters, 1.0e-8);

      // compute the residual
      Epetra_Vector bcomp( jacobian->OperatorRangeMap() );
      TEUCHOS_ASSERT_EQUALITY( 0, jacobian->Apply(*epetra_x, bcomp) );

      Epetra_Vector resid( jacobian->OperatorRangeMap() );
      TEUCHOS_ASSERT_EQUALITY( 0, resid.Update(1.0, *epetra_b, -1.0, bcomp, 0.0 ) );

      double residual;
      TEUCHOS_ASSERT_EQUALITY( 0, resid.Norm2(&residual) );
      if (eComm->MyPID()==0) cout << "Residual    = " << residual << endl;
      // -----------------------------------------------------------------------
    }
    catch ( std::exception & e )
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << e.what() << std::endl;
        status += 10;
    }
    catch ( std::string & e )
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << e << std::endl;
        status += 10;
    }
    catch ( const char * e )
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << e << std::endl;
        status += 10;
    }
    catch ( int e )
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << "Caught unknown exception code " << e <<  "." << std::endl;
        status += 10;
    }
    catch (...)
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << "Caught unknown exception." << std::endl;
        status += 10;
    }

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return status==0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
