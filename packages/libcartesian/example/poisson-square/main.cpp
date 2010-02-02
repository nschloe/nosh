// #include <Teuchos_DefaultComm.hpp>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <Thyra_LinearOpBase.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_Vector.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>

#include "GridUniformSquare.h"


// =============================================================================
// Helper function to compute a single norm for a vector
double epetraNorm2 ( const Epetra_Vector &v )
{
    double norm[1] = { -1.0 };
    v.Norm2 ( &norm[0] );
    return norm[0];
}
// =============================================================================
int main ( int argc, char *argv[] )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( &argc,&argv );
#endif

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm
    = Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm>  eComm
    = Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif


    // ===========================================================================
    // handle command line arguments
    Teuchos::CommandLineProcessor cmdProcessor;

    cmdProcessor.setDocString (
        "This program solves the Poisson problem on the square with homogenous boundary conditions.\n"
    );

    std::string extraParamsFile = "";
    cmdProcessor.setOption ( "extra", &extraParamsFile, "XML file with extra solver parameters", false );
    cmdProcessor.recogniseAllOptions ( true );  // print warning for unrecognized arguments
    cmdProcessor.throwExceptions ( true );      // do throw exceptions
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;     // finally, parse the command line
    try
    {
        parseReturn = cmdProcessor.parse ( argc, argv );
    }
    catch ( std::exception &e )
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    if ( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED )
    {
        return 0;
    }
    if ( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL )
    {
        return 1; // Error!
    }
    // =========================================================================

    // create the domain
    unsigned int nx = 100;
    double       scaling = 1.0;
    GridUniformSquare Omega ( nx, scaling );

    double h = Omega.getUniformH();

    // create Map for the domain
    int N = Omega.getNumGridPoints();
    const Epetra_Map Map ( N, 0, *eComm );

    // create matrix with the map
    // set it to be the identity matrix for now
    int numEntriesPerRow = 0;
    Teuchos::RCP<Epetra_CrsMatrix> AEpetra = Teuchos::rcp ( new Epetra_CrsMatrix ( Copy, Map, numEntriesPerRow ) );

    // construct right hand side
    Teuchos::RCP<Epetra_Vector> bEpetra = Teuchos::rcp ( new Epetra_Vector ( Map ) );


    std::vector<int>    columnsBoundary ( 1 );
    std::vector<double> valuesBoundary ( 1 );
    valuesBoundary[0] = 1.0;

    std::vector<int>    columnsInterior ( 5 );
    std::vector<double> valuesInterior ( 5 );
    valuesInterior[0] =  4.0 / ( h*h ) ;
    valuesInterior[1] = -1.0 / ( h*h ) ;
    valuesInterior[2] = -1.0 / ( h*h ) ;
    valuesInterior[3] = -1.0 / ( h*h ) ;
    valuesInterior[4] = -1.0 / ( h*h ) ;

    double rhsInterior = 1.0;
    double rhsBoundary = 0.0;
    for ( int row=0; row<Map.NumMyElements(); row++ )
    {
        int k = Map.GID ( row );
        switch ( Omega.getNodeType ( k ) )
        {
        case GridVirtual::INTERIOR:
            columnsInterior[0] = k;
            columnsInterior[1] = Omega.getKLeft ( k );
            columnsInterior[2] = Omega.getKRight ( k );
            columnsInterior[3] = Omega.getKAbove ( k );
            columnsInterior[4] = Omega.getKBelow ( k );

            AEpetra->InsertMyValues ( row, columnsInterior.size(), &valuesInterior[0], &columnsInterior[0] );
            bEpetra->ReplaceMyValue ( row, 0, rhsInterior );
            break;
            // boundary nodes
        case GridVirtual::BOTTOMLEFTCONVEX:
        case GridVirtual::BOTTOMLEFTCONCAVE:
        case GridVirtual::BOTTOMRIGHTCONVEX:
        case GridVirtual::BOTTOMRIGHTCONCAVE:
        case GridVirtual::TOPLEFTCONVEX:
        case GridVirtual::TOPLEFTCONCAVE:
        case GridVirtual::TOPRIGHTCONVEX:
        case GridVirtual::TOPRIGHTCONCAVE:
        case GridVirtual::BOTTOM:
        case GridVirtual::RIGHT:
        case GridVirtual::TOP:
        case GridVirtual::LEFT:
            columnsBoundary[0] = k;
            AEpetra->InsertMyValues ( row, columnsBoundary.size(), &valuesBoundary[0], &columnsBoundary[0] );
            bEpetra->ReplaceMyValue ( row, 0, rhsBoundary );
            break;
        }

    }
    // finalize the matrix
    AEpetra->FillComplete();


    // solve the equation system using Thyra + Stratimikos
    Teuchos::RCP<Epetra_Vector> xEpetra = Teuchos::rcp ( new Epetra_Vector ( Map ) );


    // ------------------------------------------------------------------------
    // wrap the Epetra objects into Thyra objects
    Teuchos::RCP<const Thyra::LinearOpBase<double> >
    A = Thyra::epetraLinearOp ( AEpetra );
    Teuchos::RCP<Thyra::VectorBase<double> >
    x = Thyra::create_Vector ( xEpetra, A->domain() );
    Teuchos::RCP<const Thyra::VectorBase<double> >
    b = Thyra::create_Vector ( bEpetra, A->range() );

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    // Reading in the solver parameters from the parameters file and/or from
    // the command line.  This was setup by the command-line options
    // set by the setupCLP(...) function above.
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    linearSolverBuilder.readParameters ( out.get() );

    // Augment parameters if appropriate
    if ( !extraParamsFile.empty() )
        Teuchos::updateParametersFromXmlFile ( "./"+extraParamsFile, &*linearSolverBuilder.getNonconstParameterList() );

    // Create a linear solver factory given information read from the
    // parameter list.
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = linearSolverBuilder.createLinearSolveStrategy ( "" );

    // Setup output stream and the verbosity level
    lowsFactory->setOStream ( out );
    lowsFactory->setVerbLevel ( Teuchos::VERB_LOW );

    // Create a linear solver based on the forward operator A
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> >
    lows = Thyra::linearOpWithSolve ( *lowsFactory, A );

    // Solve the linear system (note: the initial guess in 'x' is critical)
    Thyra::SolveStatus<double>
    status = Thyra::solve ( *lows, Thyra::NOTRANS, *b, &*x );
    *out << "\nSolve status:\n" << status;
    // ------------------------------------------------------------------------
    
    
    // Wipe out the Thyra wrapper for x to guarantee that the solution will be
    // written back to epetra_x!  At the time of this writting this is not
    // really needed but the behavior may change at some point so this is a
    // good idea.
    x = Teuchos::null;

    *out << "\nSolution ||epetra_x||2 = " << epetraNorm2 ( *xEpetra ) << "\n";

    *out << "\nTesting the solution error ||b-A*x||/||b|| computed through the Epetra objects ...\n";

    // r = b - A*x
    Epetra_Vector rEpetra ( *bEpetra );
    {
        Epetra_Vector AxEpetra ( AEpetra->OperatorRangeMap() );
        AEpetra->Apply ( *xEpetra,AxEpetra );
        rEpetra.Update ( -1.0,AxEpetra,1.0 );
    }

    const double tol = 1.0e-10;
    const double nrm_r = epetraNorm2 ( rEpetra );
    const double nrm_b = epetraNorm2 ( *bEpetra );
    const double rel_err = ( nrm_r / nrm_b );
    const bool   passed = ( rel_err <= tol );

    *out
    << "||b-A*x||/||b|| = " << nrm_r << "/" << nrm_b << " = " << rel_err
    << " < tol = " << tol << " ? " << ( passed ? "passed" : "failed" ) << "\n";
    
    // plot the result on the grid
    Teuchos::ParameterList params;
    std::string filePath = "test.vtk";
    Omega.writeWithGrid ( *xEpetra, params, filePath );


#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
// =========================================================================
