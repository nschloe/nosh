#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include <boost/filesystem.hpp>

#include "Ginla_EpetraFVM_ModelEvaluator.h"
#include "Ginla_EpetraFVM_State.h"
#include "VIO_EpetraMesh_Reader.h"
#include "Ginla_EpetraFVM_KineticEnergyOperator.h"

#include "Ginla_MagneticVectorPotential_X.h"
#include "Ginla_MagneticVectorPotential_Y.h"
#include "Ginla_MagneticVectorPotential_Z.h"

//#include "Ginla_IO_StateWriter.h"
//#include "Ginla_IO_StatsWriter.h"

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
      Teuchos::RCP<VIO::EpetraMesh::Mesh> mesh = Teuchos::null;

      VIO::EpetraMesh::read( eComm,
                             inputFileName,
                             z,
                             mesh,
                             problemParameters
                           );

      double mu = 1.0;
      Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp =
              Teuchos::rcp ( new Ginla::MagneticVectorPotential::Z ( mu ) );

      // create the kinetic energy operator
      Teuchos::RCP<Ginla::EpetraFVM::KineticEnergyOperator> keo =
              Teuchos::rcp( new  Ginla::EpetraFVM::KineticEnergyOperator( mesh,  mvp ) );
      // -----------------------------------------------------------------------
      // create initial guess and right-hand side
      Epetra_Vector x( keo->OperatorDomainMap() );
      Epetra_Vector b( keo->OperatorRangeMap() );

      // build the AztecOO problem
      Epetra_LinearProblem problem( &*keo, &x, &b );
      AztecOO solver( problem );
      // make sure the problem is symmetric
      problem.AssertSymmetric();

      //solver.SetAztecOption(AZ_precond, AZ_none);
      //solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
      //solver.SetAztecOption(AZ_solver, AZ_gmres);
      solver.SetAztecOption(AZ_solver, AZ_cg);
      //solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
      solver.SetAztecOption(AZ_precond, AZ_none);
      //solver.SetAztecOption(AZ_scaling, 8);
      ///solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
      //solver.SetAztecOption(AZ_output, 1);
      //solver.SetAztecOption(AZ_reorder, 0);
      solver.SetAztecOption(AZ_graph_fill, 3);
      solver.SetAztecOption(AZ_overlap, 0);
      //solver.SetAztecOption(AZ_poly_ord, 9);
      solver.SetAztecParam(AZ_ilut_fill, 4.0);
      solver.SetAztecParam(AZ_drop, 0.0);
      //double rthresh = 1.4;
      //cout << "Rel threshold = " << rthresh << endl;
      //solver.SetAztecParam(AZ_rthresh, rthresh);
      //double athresh = 10.0;
      //cout << "Abs threshold = " << athresh << endl;
      //solver.SetAztecParam(AZ_athresh, athresh);
      solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);

      int Niters = 320;
      solver.SetAztecOption(AZ_kspace, Niters);

      // do the iteration
      solver.Iterate(Niters, 1.0e-14);

      // compute the residual
      Epetra_Vector bcomp( keo->OperatorRangeMap() );
      TEUCHOS_ASSERT_EQUALITY( 0, keo->Apply(x, bcomp) );

      Epetra_Vector resid( keo->OperatorRangeMap() );
      TEUCHOS_ASSERT_EQUALITY( 0, resid.Update(1.0, b, -1.0, bcomp, 0.0 ) );

      double residual;
      TEUCHOS_ASSERT_EQUALITY( 0, resid.Norm2(&residual) );
      if (eComm->MyPID()==0) cout << "Residual    = " << residual << endl;

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
