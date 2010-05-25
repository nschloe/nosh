#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// #include <Piro_Epetra_NOXSolver.hpp>
#include <Piro_Epetra_LOCASolver.hpp>
// 
#include <boost/filesystem.hpp>

// #include "Ginla_IO_SaveNewtonData.h"
#include "Ginla_IO_NoxObserver.h"
#include "Ginla_ModelEvaluator_Default.h"
#include "Ginla_ModelEvaluator_Bordered.h"
#include "Recti_Grid_Uniform.h"
#include "Recti_Grid_Reader.h"
#include "Ginla_Komplex_LinearProblem.h"
#include "Ginla_Operator_BCCentral.h"
// #include "Ginla_Operator_BCInner.h"
// #include "Ginla_Operator_BCOuter.h"
#include "Ginla_MagneticVectorPotential_Centered.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"
#include "Ginla_IO_SaveEigenData.h"

// =============================================================================
int main ( int argc, char *argv[] )
{ 
  // Initialize MPI and timer
  int Proc=0;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  double total_time = -MPI_Wtime();
  (void) MPI_Comm_rank(MPI_COMM_WORLD, &Proc);
  MPI_Comm appComm = MPI_COMM_WORLD;
#else
  int appComm=0;
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
        "This program solves the Ginzburg--Landau problem with a Piro LOCA interface.\n"
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
    Teuchos::ParameterList outputList = paramList->sublist ( "Output", true );
    
    // set default directory to be the directory of the XML file itself
    std::string xmlPath = boost::filesystem::path ( xmlInputFileName ).branch_path().string();
    boost::filesystem::path outputDirectory = outputList.get<string> ( "Output directory" );
    if ( outputDirectory.root_directory().empty() ) // outputDirectory is empty or is a relative directory.
        outputDirectory = xmlPath / outputDirectory;
    std::string contFileBaseName =
        outputList.get<string> ( "Continuation file base name" );
    std::string outputFormat =
        outputList.get<string> ( "Output format" );
    boost::filesystem::path contDataFile =
        outputList.get<string> ( "Continuation data file name" );
        
    Teuchos::ParameterList initialGuessList;
    initialGuessList = paramList->sublist ( "Initial guess", true );
    boost::filesystem::path inputGuessFile = initialGuessList.get<string> ( "File name" );
    if ( !inputGuessFile.empty() && inputGuessFile.root_directory().empty() ) // if inputGuessFile is a relative path
        inputGuessFile = xmlPath / inputGuessFile;
    TEUCHOS_ASSERT( !inputGuessFile.empty() );
    // =========================================================================

    Teuchos::ParameterList             problemParameters;
    Teuchos::RCP<Ginla::State>         state;
    Teuchos::RCP<Recti::Grid::Uniform> grid = Teuchos::null;

    Recti::Grid::Reader::read ( Comm,
                                inputGuessFile.string(),
                                state,
                                grid,
                                problemParameters );

    // possibly overwrite the parameters
    Teuchos::ParameterList & overwriteParamsList = paramList->sublist ( "Overwrite parameter list", true ); 
    bool overwriteParameters = overwriteParamsList.get<bool> ( "Overwrite parameters" );
    if ( overwriteParameters )
    {
        Teuchos::ParameterList & overwritePList = overwriteParamsList.sublist( "Parameters", true );
        problemParameters.setParameters( overwritePList );
        
        // possibly update the scaling of the grid
        grid->updateScaling( problemParameters.get<double>("scaling") );
    }
    
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( problemParameters.get<double> ( "H0" ),
                                                                      problemParameters.get<double> ( "scaling" ) ) );

    Teuchos::RCP<Ginla::Komplex::LinearProblem> komplex =
        Teuchos::rcp( new Ginla::Komplex::LinearProblem( eComm, state->getPsi()->getMap() ) );
        
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
                                                        state->getPsi()->getMap(),
                                                        state->getPsi()->getMap() ) );

    // create the mode evaluator
//     Teuchos::RCP<Ginla::ModelEvaluator::Default> glModel = 
//               Teuchos::rcp(new Ginla::ModelEvaluator::Default( glOperator,
//                                                                komplex,
//                                                                *state,
//                                                                problemParameters ) );
    Teuchos::RCP<Ginla::ModelEvaluator::Bordered> glModel = 
              Teuchos::rcp(new Ginla::ModelEvaluator::Bordered( glOperator,
                                                                komplex,
                                                                *state,
                                                                problemParameters ) );
                                                               
    Teuchos::RCP<Ginla::IO::NoxObserver> observer =
        Teuchos::rcp( new Ginla::IO::NoxObserver( stateWriter,
                                                  glModel,
                                                  Ginla::IO::NoxObserver::CONTINUATION ) );
                                                  
    std::string contFilePath = (outputDirectory / contDataFile).string();
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter =
        Teuchos::rcp( new Ginla::IO::StatsWriter( contFilePath ) );
    observer->setStatisticsWriter( statsWriter, glOperator );

    Teuchos::RCP<Teuchos::ParameterList> piroParams =
        Teuchos::rcp(new Teuchos::ParameterList("Piro Parameters"));
    Teuchos::updateParametersFromXmlFile(xmlInputFileName, piroParams.get());

    // Use these two objects to construct a Piro solved application 
    //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
    Teuchos::RCP<EpetraExt::ModelEvaluator> piro;

//     std::string& solver = piroParams->get( "Piro Solver", "" );

//     if (solver=="NOX")
//     piro = Teuchos::rcp(new Piro::Epetra::NOXSolver( piroParams,
//                                                       glModel,
//                                                       observer ));


    // handle eigendata stuff
#ifdef HAVE_LOCA_ANASAZI
    Teuchos::RCP<Teuchos::ParameterList> eigenList = Teuchos::rcpFromRef ( piroParams->sublist ( "LOCA" ).sublist ( "Stepper" ) .sublist ( "Eigensolver" ) );
    std::string eigenvaluesFileName =
        outputDirectory.string()  + "/" + outputList.get<string> ( "Eigenvalues file name" );
    std::string eigenstateFileNameAppendix =
        outputList.get<string> ( "Eigenstate file name appendix" );

    Teuchos::RCP<Ginla::IO::StatsWriter> eigenStatsWriter =
        Teuchos::rcp( new Ginla::IO::StatsWriter( eigenvaluesFileName ) );

    Teuchos::RCP<Ginla::IO::SaveEigenData> glEigenSaver =    
        Teuchos::RCP<Ginla::IO::SaveEigenData> ( new Ginla::IO::SaveEigenData ( eigenList,
                                                                                glModel,
                                                                                stateWriter,
                                                                                eigenStatsWriter ) );

    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy =
        glEigenSaver;
    eigenList->set ( "Save Eigen Data Method", "User-Defined" );
    eigenList->set ( "User-Defined Save Eigen Data Name", "glSaveEigenDataStrategy" );
    eigenList->set ( "glSaveEigenDataStrategy", glSaveEigenDataStrategy );
#endif

      // fetch the stepper
      Teuchos::RCP<Piro::Epetra::LOCASolver> piroLOCASolver = 
          Teuchos::rcp(new Piro::Epetra::LOCASolver( piroParams,
                                                     glModel,
                                                     observer ));

      // get stepper and inject it into the eigensaver
      Teuchos::RCP<LOCA::Stepper> stepper = piroLOCASolver->getLOCAStepperNonConst();
#ifdef HAVE_LOCA_ANASAZI
        glEigenSaver->setLocaStepper ( stepper );
#endif

      piro = piroLOCASolver;
                                                      
      // Now the (somewhat cumbersome) setting of inputs and outputs
      EpetraExt::ModelEvaluator::InArgs inArgs = piro->createInArgs();
      int num_p = inArgs.Np();     // Number of *vectors* of parameters
      Teuchos::RCP<Epetra_Vector> p1 =
          Teuchos::rcp(new Epetra_Vector(*(piro->get_p_init(0))));
      inArgs.set_p(0,p1);   
      
//       RCP<Epetra_Vector> p1 = rcp(new Epetra_Vector(*(piro->get_p_init(0))));
//       int numParams = p1->MyLength(); // Number of parameters in p1 vector
//       inArgs.set_p(0,p1);


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



//       // Set output arguments to evalModel call
//       EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();
//       int num_g = outArgs.Ng(); // Number of *vectors* of responses
//       Teuchos::RCP<Epetra_Vector> g1 = Teuchos::rcp(new Epetra_Vector(*(piro->get_g_map(0))));
//       outArgs.set_g(0,g1);
//       // Solution vector is returned as extra respons vector
//       Teuchos::RCP<Epetra_Vector> gx = Teuchos::rcp(new Epetra_Vector(*(piro->get_g_map(1))));
//       outArgs.set_g(1,gx);
// 
//       Teuchos::RCP<Epetra_MultiVector> dgdp = Teuchos::rcp(new Epetra_MultiVector(g1->Map(), numParams));
//       if (computeSens) outArgs.set_DgDp(0, 0, dgdp);

      // Now, solve the problem and return the responses
      piro->evalModel(inArgs, outArgs);
      
      // manually release LOCA stepper
#ifdef HAVE_LOCA_ANASAZI
        glEigenSaver->releaseLocaStepper ();
#endif
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
