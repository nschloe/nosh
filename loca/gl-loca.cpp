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

	bool reverse = false;
	My_CLP.setOption("reverse", "forward", &reverse,
			"Orientation of the continuation in the first step");

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
	Teuchos::RCP<ComplexVector> psiLexicographic;

	if (withInitialGuess) {
		Teuchos::RCP<Tpetra::MultiVector<double, int> > psiLexicographicSplit;

		Teuchos::RCP<IoVirtual> fileIo = Teuchos::RCP<IoVirtual>(
				IoFactory::createFileIo(inputGuessFile));
		try {
			fileIo->read(Comm, psiLexicographicSplit, glParameters);
		} catch (const std::exception &e) {
			std::cout << e.what() << std::endl;
			return 1;
		}

		psiLexicographic = Teuchos::rcp(new ComplexVector(
				psiLexicographicSplit->getMap()));

		// TODO Replace this with get1dView or get2dView of the full MultiVector
		DoubleArrayRCP psiLexicographicSplitAbsView =
				psiLexicographicSplit->getVector(0)->get1dView();
		DoubleArrayRCP psiLexicographicSplitArgView =
				psiLexicographicSplit->getVector(1)->get1dView();

		int localLength = psiLexicographic->getLocalLength();
		for (int k = 0; k < localLength; k++) {
			int kGlobal = psiLexicographic->getMap()->getGlobalElement(k);
			double_complex value = polar(sqrt(
					psiLexicographicSplitAbsView[kGlobal]),
					psiLexicographicSplitArgView[kGlobal]);
			psiLexicographic->replaceLocalValue(k, value);
		}
	} else {
		Teuchos::ParameterList& glList = paramList->sublist("GL", true);
		int Nx = glList.get<int> ("Nx");
		double edgelength = glList.get<double> ("edge length");
		double H0 = glList.get<double> ("H0");

		glParameters.set("Nx", Nx);
		glParameters.set("edge length", edgelength);
		glParameters.set("H0", H0);
	}
	// ---------------------------------------------------------------------------

	// create the gl problem
	Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions =
			Teuchos::rcp(new GlBoundaryConditionsCentral());

	Teuchos::RCP<Grid> grid = Teuchos::rcp(new Grid(
			glParameters.get<int> ("Nx"), glParameters.get<double> (
					"edge length")));

	// true if no initial guess was given
	if (psiLexicographic.is_null()) {
		int NumGlobalUnknowns = grid->getNumGridPoints();
		Teuchos::RCP<Tpetra::Map<int> > standardMap = Teuchos::rcp(
				new Tpetra::Map<int>(NumGlobalUnknowns, 0, Comm));
		psiLexicographic = Teuchos::rcp(new ComplexVector(standardMap));
	}

	Teuchos::RCP<MagneticVectorPotential> A = Teuchos::rcp(
			new MagneticVectorPotential(glParameters.get<double> ("H0"),
					glParameters.get<double> ("edge length")));

	GinzburgLandau glProblem = GinzburgLandau(grid, A, boundaryConditions);

	// ---------------------------------------------------------------------------
	Teuchos::RCP<GlSystem> glsystem;
	if (withInitialGuess) {
		// If there was is an initial guess, make sure to get the ordering correct.
		// TODO:
		// Look into having this done by Trilinos. If executed on a multiproc
		// environment, we don't want p to be fully present on all processors.
		int NumComplexUnknowns = glProblem.getGrid()->getNumGridPoints();
		std::vector<int> p(NumComplexUnknowns);
		// fill p:
		glProblem.getGrid()->lexicographic2grid(&p);
		Teuchos::RCP<ComplexVector> psi = Teuchos::rcp(new ComplexVector(
				psiLexicographic->getMap()));
		// TODO:
		// The following is certainly not multiproc.
		Teuchos::ArrayRCP<const double_complex> psiView =
				psiLexicographic->get1dView();
		for (int k = 0; k < NumComplexUnknowns; k++)
			psi->replaceGlobalValue(p[k], psiView[k]);

		// Create the interface between NOX and the application
		// This object is derived from NOX::Epetra::Interface
		glsystem = Teuchos::rcp(new GlSystem(glProblem, eComm, psi,
				outputDirectory, contFileBaseName, contFileFormat,
				contDataFileName));
	} else
		glsystem = Teuchos::rcp(new GlSystem(glProblem, eComm, outputDirectory,
				contFileBaseName, contFileFormat, contDataFileName));
	// ---------------------------------------------------------------------------


	// ---------------------------------------------------------------------------
	// Create the necessary objects
	// ---------------------------------------------------------------------------
	// Create Epetra factory
	Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory = Teuchos::rcp(
			new LOCA::Epetra::Factory);

	// Create global data object
	Teuchos::RCP<LOCA::GlobalData> globalData = LOCA::createGlobalData(
			paramList, epetraFactory);

	// set the directory to which all output gets written
	glsystem->setOutputDir(outputDirectory);

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
	locaParams.addParameter("edge length", glParameters.get<double> (
			"edge length"));

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
	Teuchos::RCP<NOX::StatusTest::NormF> normF = Teuchos::rcp(
			new NOX::StatusTest::NormF(1.0e-13));
	Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters = Teuchos::rcp(
			new NOX::StatusTest::MaxIters(100));
	Teuchos::RCP<NOX::StatusTest::Generic> comboOR = Teuchos::rcp(
			new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, normF,
					maxIters));
	// ---------------------------------------------------------------------------


	Teuchos::ParameterList& bifList = paramList->sublist("LOCA").sublist(
			"Bifurcation");
	Teuchos::RCP<NOX::Abstract::Vector> lengthNormVec = Teuchos::rcp(
			new NOX::Epetra::Vector(*noxSoln));
	lengthNormVec->init(1.0);
	bifList.set("Length Normalization Vector", lengthNormVec);

	Teuchos::RCP<NOX::Abstract::Vector> initialNullVec = Teuchos::rcp(
			new NOX::Epetra::Vector(*noxSoln));
	initialNullVec->init(1.0);
	bifList.set("Initial Null Vector", initialNullVec);
	//  Teuchos::writeParameterListToXmlFile( *paramList, "input.xml");

	// ---------------------------------------------------------------------------
	// Create the stepper
	LOCA::Stepper stepper(globalData, grp, comboOR, paramList);
	// ---------------------------------------------------------------------------


	// ---------------------------------------------------------------------------
	// Perform continuation run
	LOCA::Abstract::Iterator::IteratorStatus status;
	try {
		status = stepper.run();
	} catch (char const* e) {
		std::cerr << "Exception raised: " << e << std::endl;
	}
	// ---------------------------------------------------------------------------

#ifdef HAVE_MPI
	MPI_Finalize();
#endif

	// Final return value (0 = successful, non-zero = failure)
	return status;
}
