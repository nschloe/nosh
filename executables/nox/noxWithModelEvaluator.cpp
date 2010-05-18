#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_ModelEvaluatorInterface.H>

#include <NOX_Epetra_MatrixFree.H>
#include <NOX_Epetra_FiniteDifference.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Solver_Factory.H>

#include <NOX_StatusTest_Combo.H>
#include <NOX_StatusTest_NormF.H>
#include <NOX_StatusTest_NormUpdate.H>
#include <NOX_StatusTest_NormWRMS.H>
#include <NOX_StatusTest_MaxIters.H>
#include <NOX_StatusTest_FiniteValue.H>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <EpetraExt_RowMatrixOut.h>

#include <Piro_Epetra_NOXSolver.hpp>

#include <boost/filesystem.hpp>

// #include "Ginla_IO_SaveNewtonData.h"
#include "Ginla_IO_NoxObserver.h"
#include "Ginla_ModelEvaluator_Default.h"
#include "Recti_Grid_Reader.h"
#include "Ginla_MagneticVectorPotential_Centered.h"
#include "Ginla_Operator_BCCentral.h"
#include "Ginla_Operator_BCInner.h"
#include "Ginla_Operator_BCOuter.h"
#include "Ginla_IO_StateWriter.h"

// =============================================================================
int main ( int argc, char *argv[] )
{ 
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( &argc,&argv );
#endif

    // Create a communicator for Tpetra objects
    const Teuchos::RCP<const Teuchos::Comm<int> > Comm =
          Teuchos::DefaultComm<int>::getComm();

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
        "This program solves the Ginzburg--Landau problem with a NOX interface.\n"
    );

    std::string xmlInputFileName = "";
    My_CLP.setOption ( "xml-input-file", &xmlInputFileName,
                       "XML file containing the parameter list", true );

    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions ( true );

    // do throw exceptions
    My_CLP.throwExceptions ( true );

    // finally, parse the command line
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
    My_CLP.parse ( argc, argv );

    Teuchos::RCP<Teuchos::ParameterList> paramList =
        Teuchos::rcp ( new Teuchos::ParameterList );
    std::cout << "Reading parameter list from \"" << xmlInputFileName << "\"."
              << std::endl;

    Teuchos::updateParametersFromXmlFile ( xmlInputFileName, paramList.get() );

    // =========================================================================
    // extract data of the parameter list
    Teuchos::ParameterList& ioList = paramList->sublist ( "IO", true );

    boost::filesystem::path xmlPath =
        boost::filesystem::path ( xmlInputFileName ).branch_path();

    boost::filesystem::path inputGuessFile = ioList.get<string> ( "Input guess" );
    if ( !inputGuessFile.empty() && inputGuessFile.root_directory().empty() ) // if inputGuessFile is not empty and is a relative path
        inputGuessFile = xmlPath / inputGuessFile;

    const std::string outputFormat = ioList.get<string> ( "Output format" );
    
    // set default directory to be the directory of the XML file itself
    boost::filesystem::path outputDirectory = ioList.get<string> ( "Output directory" );
    if ( outputDirectory.root_directory().empty() )
        // outputDirectory is empty or is a relative directory.
        outputDirectory = xmlPath / outputDirectory;

//     bool computeEigenvalues = paramList->sublist ( "Eigenvalues",true )
//                               .get<bool>( "Compute Eigenvalues" );
// 
//     bool computeConditionNumbers = paramList->sublist ( "Condition Numbers",true )
//                                    .get<bool> ( "Compute Condition Numbers" );
// 
//     bool plotEachNewtonStep = paramList->sublist ( "IO",true )
//                               .get<bool> ( "Plot each Newton step" );

//     std::string jacFilename = paramList->sublist ( "IO",true )
//                               .get<std::string> ( "Jacobian MATLAB matrix file name" );
    // =========================================================================

    Teuchos::ParameterList             problemParameters;
    Teuchos::RCP<Ginla::State>         state;
    Teuchos::RCP<Recti::Grid::Uniform> grid = Teuchos::null;

    Recti::Grid::Reader::read ( Comm,
                                inputGuessFile.string(),
                                state,
                                grid,
                                problemParameters );

    double h0      = 0.8; //problemParameters.get<double> ( "H0" );
    double scaling = 7.0; //problemParameters.get<double> ( "scaling" );
    
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( h0,
                                                                      scaling ) );

    Teuchos::RCP<Ginla::Komplex> komplex =
        Teuchos::rcp( new Ginla::Komplex( eComm, state->getPsiNonConst()->getMap() ) );
        
    // setup the data output
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter =
        Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory.string(),
                                                  "solution",
                                                  "VTI",
                                                  1000 ) );
                                                                      
    // create the operator
    Teuchos::RCP<Ginla::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::Operator::BCCentral ( grid,
                                                        A,
                                                        state->getPsiNonConst()->getMap(),
                                                        state->getPsiNonConst()->getMap() ) );

    // Create the interface between NOX and the application
    // This object is derived from NOX::Epetra::Interface
    Teuchos::RCP<Ginla::ModelEvaluator::Default> glModel = 
              Teuchos::rcp(new Ginla::ModelEvaluator::Default( glOperator,
                                                               komplex ) );

    Teuchos::RCP<Ginla::IO::NoxObserver> observer =
        Teuchos::rcp( new Ginla::IO::NoxObserver( stateWriter,
                                                  glModel ) );
                                                               
    Teuchos::RCP<Teuchos::ParameterList> piroParams =
        Teuchos::rcp(new Teuchos::ParameterList("Piro Parameters"));
    Teuchos::updateParametersFromXmlFile(xmlInputFileName, piroParams.get());

    // Use these two objects to construct a Piro solved application 
    //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
    Teuchos::RCP<EpetraExt::ModelEvaluator> piro;

    std::string & solver = piroParams->get( "Piro Solver", "" );
    
    std::cout << "Solver " << solver << std::endl;

//     if (solver=="NOX")
    piro = Teuchos::rcp(new Piro::Epetra::NOXSolver( piroParams,
                                                     glModel,
                                                     observer ));
      
//       bool computeSens = piroParams->get("Compute Sensitivities", false);

      // Now the (somewhat cumbersome) setting of inputs and outputs
      EpetraExt::ModelEvaluator::InArgs inArgs = piro->createInArgs();
      int num_p = inArgs.Np();     // Number of *vectors* of parameters

      Teuchos::RCP<Epetra_Vector> p1 =
          Teuchos::rcp(new Epetra_Vector(*(piro->get_p_init(0))));
      int numParams = p1->MyLength(); // Number of parameters in p1 vector
      inArgs.set_p(0,p1);

      // Set output arguments to evalModel call
      EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();
//       int num_g = outArgs.Ng(); // Number of *vectors* of responses
//       Teuchos::RCP<Epetra_Vector> g1 =
//           Teuchos::rcp(new Epetra_Vector(*(piro->get_g_map(0))));
//       outArgs.set_g(0,g1);
      // Solution vector is returned as extra respons vector
//       Teuchos::RCP<Epetra_Vector> gx =
//           Teuchos::rcp(new Epetra_Vector(*(piro->get_g_map(1))));
//       outArgs.set_g(1,gx);

//       Teuchos::RCP<Epetra_MultiVector> dgdp =
//           Teuchos::rcp(new Epetra_MultiVector(g1->Map(), numParams));
//       if (computeSens) outArgs.set_DgDp(0, 0, dgdp);

      // Now, solve the problem and return the responses
      piro->evalModel(inArgs, outArgs);                                    

    }
    catch ( std::exception & e )
    {
        std::cerr << e.what() << std::endl;
        status += 10;
    }
    catch (...)
    {
        std::cerr << "Caught unknown exception." << std::endl;
        status += 10;
    }
      
#ifdef HAVE_MPI
      MPI_Finalize();
#endif

      // Final return value (0 = successfull, non-zero = failure)
      return status;
}
// =========================================================================
