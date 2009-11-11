#include <LOCA.H>
#include <LOCA_Epetra_Factory.H>
#include <LOCA_Epetra_Group.H>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Vector.h>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "ioFactory.h"
#include "ioException.h"

#include "glSystem.h"
#include "glBoundaryConditionsInner.h"
#include "glBoundaryConditionsOuter.h"
#include "glBoundaryConditionsCentral.h"
#include "eigenSaver.h"

typedef complex<double> double_complex;
typedef Tpetra::Vector<double_complex, int> ComplexVector;
typedef Teuchos::ArrayRCP<const double> DoubleArrayRCP;
typedef Teuchos::ArrayRCP<const double_complex> ComplexArrayRCP;

int main(int argc, char *argv[]) {

	// Initialize MPI
#ifdef HAVE_MPI
	MPI_Init(&argc,&argv);
#endif

	// create Epetra communicator
#ifdef HAVE_MPI
	Teuchos::RCP<Epetra_MpiComm> eComm =
	Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
	Teuchos::RCP<Epetra_SerialComm> eComm = Teuchos::rcp<Epetra_SerialComm>(
			new Epetra_SerialComm());
#endif

	// Create a communicator for Tpetra objects
	const Teuchos::RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<
			int>::getComm();

	// =========================================================================
	// handle command line arguments
	Teuchos::CommandLineProcessor My_CLP;

	My_CLP.setDocString(
			"This program does continuation for the Ginzburg--Landau problem with a LOCA interface.\n"
				"It is possible to give an initial guess in VTK format on the command line.\n");

//	bool reverse = false;
//	My_CLP.setOption("reverse", "forward", &reverse,
//			"Orientation of the continuation in the first step");

	//	bool stopOnUnstable = false;
	//	My_CLP.setOption("stop-on-unstable", "continue-on-unstable",
	//			&stopOnUnstable,
	//			"[unused] Stop the continuation if an unstable state has been reached");

	std::string xmlInputFileName = "";
	My_CLP.setOption("xml-input-file", &xmlInputFileName,
			"XML file containing the parameter list", true);

	// print warning for unrecognized arguments
	My_CLP.recogniseAllOptions(true);

	// don't throw exceptions
	My_CLP.throwExceptions(false);

	// finally, parse the stuff!
	Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn =
			My_CLP.parse(argc, argv);
	if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
		return 0;
	}
	if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
		return 1; // Error!
	}

	Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(
			new Teuchos::ParameterList);
	std::cout << "Reading parameter list from \"" << xmlInputFileName << "\"."
			<< std::endl;
	Teuchos::updateParametersFromXmlFile(xmlInputFileName, paramList.get());
	// =========================================================================
	// extract data of the parameter list
	Teuchos::ParameterList& ioList = paramList->sublist("IO", true);
	std::string inputGuessFile = ioList.get<string> ("Input guess");
	bool withInitialGuess = inputGuessFile.length() > 0;
	std::string outputDirectory = ioList.get<string> ("Output directory");
	std::string contFileBaseName = ioList.get<string> (
			"Continuation file base name");
	std::string contFileFormat =
			ioList.get<string> ("Continuation file format");
	std::string contDataFileName = ioList.get<string> (
			"Continuation data file name");
	// =========================================================================


	// ---------------------------------------------------------------------------
	// define a new dummy psiLexicographic vector, to be adapted instantly
	Teuchos::ParameterList glParameters;
	Teuchos::RCP<ComplexVector> psi;
	Teuchos::RCP<Grid> grid;

	if (withInitialGuess) {
	   try {
		   readStateFromFile( Comm, inputGuessFile, psi, grid, glParameters );
	   }
	   catch (const IoException& e) {
	       std::cerr << e.what() << std::endl;
	   }
	   catch (...) {
	       std::cerr << "Unknown exception caught." << std::endl;
	   }
	} else {
		// read the parameters from the XML file
		Teuchos::ParameterList& glList = paramList->sublist("GL", true);
		int Nx = glList.get<int> ("Nx");
		double edgeLength = glList.get<double> ("edge length");
		double H0 = glList.get<double> ("H0");

		glParameters.set("Nx", Nx);
		glParameters.set("edge length", edgeLength);
		glParameters.set("H0", H0);

		grid = Teuchos::rcp(new Grid(Nx,edgeLength) );

		// set initial guess
		int numComplexUnknowns = grid->getNumGridPoints();
		Teuchos::RCP<Tpetra::Map<int> >complexMap =
				   Teuchos::rcp(new Tpetra::Map<int>(numComplexUnknowns,0,Comm) );
		psi = Teuchos::rcp( new ComplexVector(complexMap) );
		double_complex alpha( 1.0, 0.0 );
		psi->putScalar( alpha );
	}
	// ---------------------------------------------------------------------------

	// Make sure that the calculation starts off with the correct H0.
	// TODO See if it's possible to declare the initial parameter once instead of three times in the
	//      input XML file.
//	paramList->sublist("LOCA").sublist("Stepper").set("Initial Value", glParameters.get<double>("H0") );

	// create the gl problem
	Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions =
			Teuchos::rcp(new GlBoundaryConditionsCentral());

	Teuchos::RCP<MagneticVectorPotential> A = Teuchos::rcp(
			new MagneticVectorPotential( glParameters.get<double> ("H0"),
				                         glParameters.get<double> ("edge length")));

	GinzburgLandau glProblem = GinzburgLandau(grid, A, boundaryConditions);

	Teuchos::RCP<GlSystem> glsystem = Teuchos::rcp(new GlSystem( glProblem, eComm, psi,
			                                                     outputDirectory, contDataFileName,
			                                                     contFileFormat, contFileBaseName ));

	// ---------------------------------------------------------------------------
	// Create the necessary objects
	// ---------------------------------------------------------------------------
	// Create Epetra factory
	Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory = Teuchos::rcp(
			new LOCA::Epetra::Factory);

	// Create global data object
	Teuchos::RCP<LOCA::GlobalData> globalData = LOCA::createGlobalData(
			paramList, epetraFactory);

	// get the initial solution
	Teuchos::RCP<Epetra_Vector> soln = glsystem->getSolution();
	// ---------------------------------------------------------------------------


#ifdef HAVE_LOCA_ANASAZI
	Teuchos::ParameterList& eigenList = paramList->sublist("LOCA").sublist(
			"Stepper") .sublist("Eigensolver");
	Teuchos::RCP<Teuchos::ParameterList> eigenListPtr = Teuchos::rcpFromRef(
			(eigenList));
	std::string eigenvaluesFileName = ioList.get<string> (
			"Eigenvalues file name");
	std::string eigenstateFileNameAppendix = ioList.get<string> (
			"Eigenstate file name appendix");
	Teuchos::RCP<EigenSaver> glEigenSaver = Teuchos::RCP<EigenSaver>(
			new EigenSaver(eigenListPtr, globalData, outputDirectory,
					eigenvaluesFileName, contFileBaseName,
					eigenstateFileNameAppendix, glsystem));

	Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>
			glSaveEigenDataStrategy = glEigenSaver;
	eigenList.set("glSaveEigenDataStrategy", glSaveEigenDataStrategy);
#endif

	// ---------------------------------------------------------------------------
	// Create all possible Epetra_Operators.
	Teuchos::RCP<Epetra_RowMatrix> Analytic = glsystem->getJacobian();

	// Create the linear system.
	// Use the TimeDependent interface for computation of shifted matrices.
	Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = glsystem;
	Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glsystem;

	Teuchos::ParameterList& nlPrintParams = paramList->sublist("NOX") .sublist(
			"Printing");

	Teuchos::ParameterList& lsParams = paramList->sublist("NOX") .sublist(
			"Direction") .sublist("Newton") .sublist("Linear Solver");

	Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(
			new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams, iReq,
					iJac, Analytic, *soln));

	Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = glsystem;
	// ---------------------------------------------------------------------------


	// ---------------------------------------------------------------------------
	// Create a group which uses that problem interface. The group will
	// be initialized to contain the default initial guess for the
	// specified problem.
	LOCA::ParameterVector locaParams;

	locaParams.addParameter("H0", glParameters.get<double> ("H0"));
	locaParams.addParameter("edge length", glParameters.get<double>("edge length") );

	NOX::Epetra::Vector initialGuess(soln, NOX::Epetra::Vector::CreateView);

	Teuchos::RCP<LOCA::Epetra::Group> grp = Teuchos::rcp(
			new LOCA::Epetra::Group(globalData, nlPrintParams, iTime,
					initialGuess, linSys, linSys, locaParams));

	grp->setParams(locaParams);
	// ---------------------------------------------------------------------------


	// ---------------------------------------------------------------------------
	// Get the vector from the Problem
	Teuchos::RCP<NOX::Epetra::Vector> noxSoln = Teuchos::rcp(
			new NOX::Epetra::Vector(soln, NOX::Epetra::Vector::CreateView));
	// ---------------------------------------------------------------------------


	// ---------------------------------------------------------------------------
	// Set up the status tests
	// ---------------------------------------------------------------------------
	Teuchos::ParameterList& noxList = paramList->sublist("NOX", true);
	double tol = noxList.get<double>("Tolerance");
	int maxNonlinarSteps = noxList.get<int>("Max steps");
	Teuchos::RCP<NOX::StatusTest::NormF> normF = Teuchos::rcp(
			new NOX::StatusTest::NormF(tol));
	Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters = Teuchos::rcp(
			new NOX::StatusTest::MaxIters(maxNonlinarSteps));
	Teuchos::RCP<NOX::StatusTest::Generic> comboOR = Teuchos::rcp(
			new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, normF,
					maxIters));
	// ---------------------------------------------------------------------------



	// ---------------------------------------------------------------------------
	// Create the stepper
	Teuchos::RCP<LOCA::Stepper> stepper = Teuchos::rcp( new LOCA::Stepper(globalData, grp, comboOR, paramList) );
	// ---------------------------------------------------------------------------

	// make sure that the stepper starts off with the correct starting value

	// pass pointer to stepper to glSystem to be able to read stats from the stepper in there
	glsystem->setLocaStepper( stepper );
#ifdef HAVE_LOCA_ANASAZI
	glEigenSaver->setLocaStepper( stepper );
#endif

	// ---------------------------------------------------------------------------
	// Perform continuation run
	LOCA::Abstract::Iterator::IteratorStatus status;
	try {
		status = stepper->run();
	} catch (char const* e) {
		std::cerr << "Exception raised: " << e << std::endl;
	}
	// ---------------------------------------------------------------------------


	LOCA::destroyGlobalData(globalData);
	stepper = Teuchos::null;
	glsystem = Teuchos::null;
#ifdef HAVE_LOCA_ANASAZI
	glEigenSaver = Teuchos::null;
#endif

#ifdef HAVE_MPI
	MPI_Finalize();
#endif

	// Final return value (0 = successful, non-zero = failure)
	return status;
}