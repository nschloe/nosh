// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Vector.h>

#include <LOCA_SaveEigenData_AbstractStrategy.H>
#include <LOCA_Abstract_Factory.H>
#include <LOCA_Epetra_Factory.H>
#include <LOCA_Epetra_ModelEvaluatorInterface.H>
#include <LOCA_Epetra_Group.H>
#include <LOCA_StatusTest_Combo.H>
#include <LOCA_StatusTest_MaxIters.H>
#include <NOX_Epetra_LinearSystem_Stratimikos.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_StatusTest_Generic.H>
#include <NOX_StatusTest_Factory.H>

#include "Ginla_StkMeshReader.hpp"

#include "Ginla_State.hpp"
#include "Ginla_ModelEvaluator.hpp"
#include "Ginla_StatsWriter.hpp"
#include "Ginla_NoxObserver.hpp"
#include "Ginla_SaveEigenData.hpp"
#include "Ginla_MagneticVectorPotential.hpp"

#include <Teuchos_TimeMonitor.hpp>

//#include<fenv.h>

// get pi
#define _USE_MATH_DEFINES

// =============================================================================
std::string
extractDirectory( const std::string& path )
{
    return path.substr( 0, path.find_last_of( '/' ) +1 );
}
// =============================================================================
int
main ( int argc, char *argv[] )
{
//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
//feenableexcept(FE_DIVBYZERO);

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
        // =====================================================================
        // handle command line arguments
        Teuchos::CommandLineProcessor My_CLP;

        My_CLP.setDocString (
                "This program solves the Ginzburg--Landau problem using LOCA.\n"
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

        Teuchos::RCP<Teuchos::ParameterList> piroParams =
                Teuchos::rcp ( new Teuchos::ParameterList );
        *out << "\n\n\tReading parameter list from \"" << xmlInputFileName << "\".\n\n"
             << std::endl;

        Teuchos::updateParametersFromXmlFile ( xmlInputFileName, piroParams.get() );

        // =====================================================================
        // extract data of the parameter list
        Teuchos::ParameterList outputList = piroParams->sublist ( "Output", true );

        // set default directory to be the directory of the XML file itself
        std::string xmlDirectory = extractDirectory( xmlInputFileName );

        std::string & outputDirectory = xmlDirectory;

        std::string contFilePath        = xmlDirectory + '/' + outputList.get<std::string> ( "Continuation data file name" );
        std::string eigenvaluesFilePath = xmlDirectory + '/' + outputList.get<std::string> ( "Eigenvalues file name" );

        Teuchos::ParameterList initialGuessList;
        initialGuessList = piroParams->sublist ( "Initial guess", true );
	std::string inputFilePath = xmlDirectory + '/' + initialGuessList.get<std::string> ( "State" );
        // =====================================================================
        // Read the data from the file.
        Teuchos::ParameterList data;
        Ginla::StkMeshRead( *eComm, inputFilePath, data );

        // Cast the data into something more accessible.
        Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
        Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
        Teuchos::RCP<const Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<const Epetra_MultiVector>() );
        Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
        Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

        // set the output directory for later plotting with this
        mesh->setOutputFile( outputDirectory, "solution" );

        // create the state
        TEUCHOS_ASSERT( !z.is_null() );
        Teuchos::RCP<Ginla::State> state =
            Teuchos::rcp( new Ginla::State( *z, mesh ) );

        // possibly overwrite the parameters
        Teuchos::ParameterList & overwriteParamsList = piroParams->sublist ( "Overwrite parameter list", true );
        bool overwriteParameters = overwriteParamsList.get<bool> ( "Overwrite parameters" );
        if ( overwriteParameters )
        {
            Teuchos::ParameterList & overwritePList = overwriteParamsList.sublist( "Parameters", true );
            problemParameters.setParameters( overwritePList );
        }

        double mu = problemParameters.get<double> ( "mu" );

        Teuchos::RCP<Ginla::MagneticVectorPotential> mvp =
                Teuchos::rcp ( new Ginla::MagneticVectorPotential ( mesh, mvpValues, mu ) );

        // create the mode evaluator
        Teuchos::RCP<Ginla::ModelEvaluator> glModel =
                Teuchos::rcp( new Ginla::ModelEvaluator( mesh,
                                                         problemParameters,
                                                         thickness,
                                                         mvp,
                                                         z
                                                       )
                            );

        // set I/O routines
        int maxLocaSteps = piroParams->sublist ( "LOCA" )
                                      .sublist ( "Stepper" )
                                      .get<int> ( "Max Steps" );

        Teuchos::RCP<Ginla::NoxObserver> observer =
                Teuchos::rcp( new Ginla::NoxObserver( glModel,
                                                      Ginla::NoxObserver::OBSERVER_TYPE_CONTINUATION
                                                    )
                            );

        Teuchos::RCP<Ginla::StatsWriter> statsWriter =
                Teuchos::rcp( new Ginla::StatsWriter( contFilePath ) );
        observer->setStatisticsWriter( statsWriter );


//         setup eigen saver
#ifdef HAVE_LOCA_ANASAZI
          Teuchos::ParameterList & eigenList = piroParams->sublist ( "LOCA" ).sublist ( "Stepper" ) .sublist ( "Eigensolver" );
//          std::string eigenstateFileNameAppendix =
//              outputList.get<std::string> ( "Eigenstate file name appendix" );
          Teuchos::RCP<Ginla::StatsWriter> eigenStatsWriter =
              Teuchos::rcp( new Ginla::StatsWriter( eigenvaluesFilePath ) );

          Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy =
                  Teuchos::RCP<Ginla::SaveEigenData> ( new Ginla::SaveEigenData ( eigenList,
                                                                                  glModel,
                                                                                  eigenStatsWriter ) );
          eigenList.set ( "Save Eigen Data Method", "User-Defined" );
          eigenList.set ( "User-Defined Save Eigen Data Name", "glSaveEigenDataStrategy" );
          eigenList.set ( "glSaveEigenDataStrategy", glSaveEigenDataStrategy );
#endif

        // explicitly copy the parameter list as LOCA::createGlobalData needs a non-const
        Teuchos::RCP<Teuchos::ParameterList> locaParams =
                Teuchos::rcp( new Teuchos::ParameterList( piroParams ->sublist ( "LOCA" ) ) );

        std::string contParam = piroParams->sublist ( "LOCA" )
                                .sublist ( "Stepper" )
                                .get<std::string> ( "Continuation Parameter" );
        if ( eComm->MyPID() == 0 && problemParameters.isParameter ( contParam ) )
            std::cerr << "Warning: Continuation parameter \""
                    << contParam
                    << "\" explicitly given. Initial value will be overwritten by 'LOCA->Stepper->Initial Value', though."
                    << std::endl;

        // ---------------------------------------------------------------------
        // Create Epetra factory
        Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
                Teuchos::rcp ( new LOCA::Epetra::Factory() );

        // Create global data object
        Teuchos::RCP<LOCA::GlobalData> globalData =
                LOCA::createGlobalData ( locaParams, epetraFactory );
        // ---------------------------------------------------------------------
        Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> locaModelEvaluatorInterface
                = Teuchos::rcp( new LOCA::Epetra::ModelEvaluatorInterface( globalData, glModel ) );
        // ---------------------------------------------------------------------
        locaModelEvaluatorInterface->setObserver( observer );
        // ---------------------------------------------------------------------
        // Create all possible Epetra_Operators.
        Teuchos::RCP<Epetra_Operator> J = glModel->create_W();
        Teuchos::RCP<Epetra_Operator> M = glModel->create_WPrec()->PrecOp;

        // Create the linear system.
        // Use the TimeDependent interface for computation of shifted matrices.
        //Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = locaModelEvaluatorInterface;
        Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = locaModelEvaluatorInterface;
        Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = locaModelEvaluatorInterface;

        Teuchos::ParameterList& nlPrintParams =
                piroParams->sublist ( "NOX" ) .sublist ( "Printing", true );

        //  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(
        //      new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams, iJac, J, iPrec, M, *soln));

        NOX::Epetra::Vector cloneVector( Epetra_Vector( *glModel->get_x_map() ) );

//        // AztecOO system
//        Teuchos::ParameterList& lsParams =  piroParams->sublist ( "NOX" )
//                                                       .sublist ( "Direction" )
//                                                       .sublist ( "Newton" )
//                                                       .sublist ( "Linear Solver", true );
//        Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
//                Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( nlPrintParams,
//                                                                      lsParams,
//                                                                      iJac,
//                                                                      J,
//                                                                      iPrec,
//                                                                      M,
//                                                                      cloneVector ) );

        // Stratimikos system
        Teuchos::ParameterList& lsParams = piroParams->sublist ( "NOX" )
                                                      .sublist ( "Direction" )
                                                      .sublist ( "Newton" )
                                                      .sublist ( "Stratimikos Linear Solver", true );
        // is this really needed?
        lsParams.set( "Preconditioner", "User Defined" );

        bool isAlreadyInverted = true;

        Teuchos::RCP<NOX::Epetra::LinearSystemStratimikos> linSys =
                Teuchos::rcp ( new NOX::Epetra::LinearSystemStratimikos ( nlPrintParams,
                                                                          lsParams,
                                                                          iJac,
                                                                          J,
                                                                          iPrec,
                                                                          M,
                                                                          cloneVector,
                                                                          isAlreadyInverted ) );

        Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = locaModelEvaluatorInterface;
        // ---------------------------------------------------------------------
        // Create a group which uses that problem interface. The group will
        // be initialized to contain the default initial guess for the
        // specified problem.
        LOCA::ParameterVector myLocaParams = *(glModel->getParameters());

        // get initial guess
        NOX::Epetra::Vector initialGuess ( state->getPsiNonConst(),
                                           NOX::Epetra::Vector::CreateView
                                           );

        Teuchos::RCP<LOCA::Epetra::Group> grp =
                Teuchos::rcp ( new LOCA::Epetra::Group ( globalData,
                                                         nlPrintParams,
                                                         iTime,
                                                         initialGuess,
                                                         linSys,
                                                         linSys,
                                                         myLocaParams ) );

        *out << myLocaParams << std::endl;
        grp->setParams ( myLocaParams );
        // ---------------------------------------------------------------------
        // Set up the NOX status tests
        Teuchos::ParameterList& noxList = piroParams->sublist ( "NOX", true );
        Teuchos::RCP<NOX::StatusTest::Generic> noxStatusTest = NOX::StatusTest::buildStatusTests( noxList.sublist( "Status Tests" ),
                                                                                                  *globalData->locaUtils()
                                                                                                  );
        // ---------------------------------------------------------------------
        // Set up the LOCA status tests
        Teuchos::RCP<LOCA::StatusTest::Combo> locaCombo =
                Teuchos::rcp ( new LOCA::StatusTest::Combo ( LOCA::StatusTest::Combo::OR ) );

        Teuchos::RCP<LOCA::StatusTest::Abstract> maxLocaStepsTest =
                Teuchos::rcp ( new LOCA::StatusTest::MaxIters ( maxLocaSteps ) );
        locaCombo->addStatusTest( maxLocaStepsTest );

        //      Teuchos::RCP<LOCA::StatusTest::Abstract> zeroEnergyTest =
        //                       Teuchos::rcp ( new Ginla::StatusTest::Energy ( glModel,
        //                                                                      0.0 ) );

        //      Teuchos::RCP<LOCA::StatusTest::Abstract> loopTest =
        //                       Teuchos::rcp ( new Ginla::StatusTest::Loop ( glsystem ) );

        //      Teuchos::RCP<LOCA::StatusTest::Abstract> turnaroundTest =
        //                       Teuchos::rcp ( new Ginla::StatusTest::Turnaround () );



        //     locaCombo->addStatusTest( zeroEnergyTest );
        //      locaCombo->addStatusTest( loopTest );
        //      locaCombo->addStatusTest( turnaroundTest );
        // ---------------------------------------------------------------------
        // Create the stepper
        Teuchos::RCP<LOCA::Stepper> stepper =
                Teuchos::rcp ( new LOCA::Stepper ( globalData,
                                                   grp,
                                                   locaCombo,
                                                   noxStatusTest,
                                                   piroParams ) );
        // ---------------------------------------------------------------------
        // pass pointer to stepper to glsystem to be able to read stats from the stepper in there
//        locaModelEvaluatorInterface->setLocaStepper ( stepper );
        // ---------------------------------------------------------------------
        // Perform continuation run
        Teuchos::RCP<Teuchos::Time> locarunTime = Teuchos::TimeMonitor::getNewTimer("LOCA runtime");
        {
            Teuchos::TimeMonitor tm(*locarunTime);
            int ierr = stepper->run();
            if ( ierr!=0 )
                success = false;
        }
        // ---------------------------------------------------------------------
        // clean up
        LOCA::destroyGlobalData ( globalData );
//        locaModelEvaluatorInterface->releaseLocaStepper();
        // ---------------------------------------------------------------------
        // print timing data
        Teuchos::TimeMonitor::summarize();
        // ---------------------------------------------------------------------
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
