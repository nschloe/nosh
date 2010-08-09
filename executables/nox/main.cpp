#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <EpetraExt_RowMatrixOut.h>

#include <boost/filesystem.hpp>

#include "glNoxHelpers.h"

#include "Ginla_IO_SaveNewtonData.h"
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
        parseReturn = My_CLP.parse ( argc, argv );

        if ( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED )
            return 0;

        if ( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL )
            return 1; // Error!

        Teuchos::RCP<Teuchos::ParameterList> paramList =
            Teuchos::rcp ( new Teuchos::ParameterList );
        if ( Comm->getRank()==0 )
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
        {
            // outputDirectory is empty or is a relative directory.
            outputDirectory = xmlPath / outputDirectory;
        }

        bool computeEigenvalues = paramList->sublist ( "Eigenvalues",true )
                                  .get<bool>( "Compute Eigenvalues" );

        bool computeConditionNumbers = paramList->sublist ( "Condition Numbers",true )
                                      .get<bool> ( "Compute Condition Numbers" );

        bool plotEachNewtonStep = paramList->sublist ( "IO",true )
                                  .get<bool> ( "Plot each Newton step" );

        std::string jacFilename = paramList->sublist ( "IO",true )
                                  .get<std::string> ( "Jacobian MATLAB matrix file name" );
        // =========================================================================

        // set problemParameters and glSystem
        Teuchos::ParameterList                    problemParameters;
        Teuchos::RCP<Ginla::FDM::LocaSystem::Virtual>  glSystem = Teuchos::null;
        Teuchos::RCP<Recti::Grid::Uniform>        grid = Teuchos::null;
        Teuchos::RCP<Ginla::FDM::State> initialState = Teuchos::null;
        if ( !inputGuessFile.empty() )
        {
            glNoxHelpers::createGlSystem ( Comm,
                                           eComm,
                                           inputGuessFile.string(),
                                           problemParameters,
                                           glSystem,
                                           initialState,
                                           grid );
        }
        else
        {
            int    Nx      = paramList->sublist ( "GL",true ).get<int> ( "Nx" );
            double scaling = paramList->sublist ( "GL",true ).get<double> ( "scaling" );
            double H0      = paramList->sublist ( "GL",true ).get<double> ( "H0" );

            paramList->sublist ( "GL",true ).get ( "chi",0.0 );
            Teuchos::ParameterList & domainParameters =
                    paramList->sublist ( "Domain",true );
            glNoxHelpers::createGlSystem ( Comm,
                                           eComm,
                                           Nx,
                                           scaling,
                                           H0,
                                           domainParameters,
                                           problemParameters,
                                           glSystem,
                                           initialState,
                                           grid );
        }
    
        Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
            Teuchos::rcpFromRef ( paramList->sublist ( "NOX",true ) );

        if ( plotEachNewtonStep )
        {
            const std::string fileBaseName = "newtonStep";
            const std::string outputFormat = "VTI";
            const unsigned int maxIndex    = 1000; // TODO replace by maxnumsteps
            Teuchos::RCP<Ginla::IO::StateWriter> stateWriter =
                Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory.string(),
                                                          fileBaseName,
                                                          outputFormat,
                                                          maxIndex ) );
                        
            Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo =
                Teuchos::rcp ( new Ginla::IO::SaveNewtonData ( stateWriter,
                                                               glSystem ) );
            nlParamsPtr->sublist ( "Solver Options" )
                        .set ( "User Defined Pre/Post Operator", ppo );
        }

        Teuchos::RCP<NOX::Epetra::Group> grpPtr =
            glNoxHelpers::createSolverGroup ( glSystem,
                                              nlParamsPtr,
                                              initialState
                                            );

        // create convergence test from sublist
        Teuchos::RCP<NOX::StatusTest::Generic> statusTest =
            glNoxHelpers::createConvergenceTest ( paramList->sublist ( "NOX Status Test",true ),
                                                  *nlParamsPtr );

        // create solver object
        Teuchos::RCP<NOX::Solver::Generic> solver =
            glNoxHelpers::createSolver ( grpPtr, statusTest, nlParamsPtr );
            
        // solve the system
        NOX::StatusTest::StatusType solvStatus = solver->solve();
        
        // compute the condition number
        if ( computeConditionNumbers )
        {
            double kappa = glNoxHelpers::computeJacobianConditionNumber ( solver, grpPtr );
            std::cout << "Condition number: kappa = " << kappa << "." << std::endl;
        }

        // spit out the Jacobian in MATLAB readable format
        if ( !jacFilename.empty() )
        {
    //       EpetraExt::RowMatrixToMatlabFile(jacFilename.c_str(),*(glsystem->getJacobian()));
        }

        // compute the eigenvalues of the Jacobian
        if ( computeEigenvalues )
            glNoxHelpers::computeJacobianEigenvalues ( solver, grpPtr, Comm->getRank() );
        
        // print the solution to a file
        const std::string fileBaseName = "solution";
        
        glNoxHelpers::printSolutionToFile ( outputDirectory.string(),
                                            fileBaseName,
                                            "VTI",
                                            solver,
                                            glSystem );
                                            
        // check the convergence status
        int status = glNoxHelpers::checkConvergence ( solver );
    }
    catch ( std::exception & e )
    {
        std::cerr << e.what() << std::endl;
        status += 10;
    }
    catch ( int & e )
    {
        std::cerr << "Caught error code \"" << e << "\"." << std::endl;
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

    return status;
}
// =========================================================================
