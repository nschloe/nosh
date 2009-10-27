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
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "ioFactory.h"
#include "glSystem.h"
#include "glBoundaryConditionsInner.h"
#include "glBoundaryConditionsOuter.h"
#include "glBoundaryConditionsCentral.h"
#include "eigenSaver.h"

typedef complex<double> double_complex;

int main(int argc, char *argv[])
{

// Initialize MPI
#ifdef HAVE_MPI
MPI_Init(&argc,&argv);
#endif

// create Epetra communicator
#ifdef HAVE_MPI
Teuchos::RCP<Epetra_MpiComm> eComm =
 Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
Teuchos::RCP<Epetra_SerialComm> eComm =
	    Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

// Create a communicator for Tpetra objects
const Teuchos::RCP<const Teuchos::Comm<int> > Comm
				 = Teuchos::DefaultComm<int>::getComm();

// =========================================================================
// handle command line arguments
Teuchos::CommandLineProcessor My_CLP;

My_CLP.setDocString(
"This program does continuation for the Ginzburg--Landau problem with a LOCA interace.\n"
"It is possible to give an initial guess in VTK format on the command line.\n"
);

bool verbose=false;
My_CLP.setOption("verbose", "silent", &verbose, "Verbostity flag" );

bool reverse=false;
My_CLP.setOption("reverse","forward",&reverse, "Orientation of the continuation in the first step" );

bool computeEigenvalues;
My_CLP.setOption("eigenvalues","no-eigenvalues",&computeEigenvalues, "Compute the eigenvalues of the Jacobian on the fly" );

std::string filename = "";
My_CLP.setOption("input-guess", &filename, "File name with initial guess");

std::string outputdir = "data";
My_CLP.setOption("output-dir", &outputdir, "Directory to which all the solution files are written");

// print warning for unrecognized arguments
My_CLP.recogniseAllOptions(true);

// don't throw exceptions
My_CLP.throwExceptions(false);

// finally, parse the stuff!
Teuchos::CommandLineProcessor::EParseCommandLineReturn
parseReturn= My_CLP.parse( argc, argv );
if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
return 0;
}
if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
return 1; // Error!
}

bool withInitialGuess = filename.length()>0;
// =========================================================================


  // ---------------------------------------------------------------------------
  Teuchos::ParameterList problemParameters;
  // define a new dummy psiLexicographic vector, to be adapted instantly
  Teuchos::RCP<Tpetra::Map<int> > dummyMap =
    Teuchos::rcp ( new Tpetra::Map<int> ( 1, 0, Comm ) );
  Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psiLexicographic =
                                                               Teuchos::ENull();

//       = Teuchos::rcp( new Tpetra::MultiVector<double_complex,int>(dummyMap,1) );
//       = Teuchos::RCP::;

  if ( withInitialGuess )
    {
      Teuchos::RCP<IoVirtual> fileIo =
               Teuchos::RCP<IoVirtual> ( IoFactory::createFileIo ( filename ) );
      try
        {
          fileIo->read ( psiLexicographic,
			 Comm,
                         problemParameters );
        }
      catch ( const std::exception &e )
        {
          std::cout << e.what() << std::endl;
          return 1;
        }
    }
  else
    {
      // set the default value
      int    Nx         = 50;
      double edgelength = 10.0;
      double H0         = 0.3;
      std::cout << "Using the standard parameters \n"
                << "    Nx         = " << Nx << ",\n"
                << "    edgelength = " << edgelength << ",\n"
                << "    H0         = " << H0 << "." << std::endl;
      problemParameters.set ( "Nx"        , Nx );
      problemParameters.set ( "edgelength", edgelength );
      problemParameters.set ( "H0"        , H0 );

int NumGlobalUnknowns = ( Nx+1 ) * ( Nx+1 );
Teuchos::RCP<Tpetra::Map<int> > standardMap
 = Teuchos::rcp ( new Tpetra::Map<int> ( NumGlobalUnknowns, 0, Comm ) );
psiLexicographic = Teuchos::rcp ( new Tpetra::MultiVector<double_complex,int> ( standardMap,1 ) );
//       psiLexicographic->replaceMap( standardMap );
}
// ---------------------------------------------------------------------------

// create the gl problem
Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions =
		       Teuchos::rcp(new GlBoundaryConditionsCentral() );

  Teuchos::RCP<StaggeredGrid> sGrid =
              Teuchos::rcp( new StaggeredGrid( problemParameters.get<int>("Nx"),
                                               problemParameters.get<double>("edgelength"),
                                               problemParameters.get<double>("H0")
                                             )
                          );

GinzburgLandau glProblem = GinzburgLandau( sGrid,
				     boundaryConditions
				   );


// ---------------------------------------------------------------------------
Teuchos::RCP<GlSystem> glsystem;
if ( withInitialGuess )
{
// If there was is an initial guess, make sure to get the ordering correct.
// TODO:
// Look into having this done by Trilinos. If executed on a multiproc
// environment, we don't want p to be fully present on all processors.
int NumComplexUnknowns = glProblem.getStaggeredGrid()->getNumComplexUnknowns();
std::vector<int> p ( NumComplexUnknowns );
// fill p:
glProblem.getStaggeredGrid()->lexicographic2grid ( &p );
Teuchos::RCP<Tpetra::MultiVector<double_complex,int> >  psi
= Teuchos::rcp ( new Tpetra::MultiVector<double_complex,int> ( psiLexicographic->getMap(),1 ) );
// TODO:
// The following is certainly not multiproc.
Teuchos::ArrayRCP<const double_complex> psiView = psiLexicographic->getVector ( 0 )->get1dView();
for ( int k=0; k<NumComplexUnknowns; k++ )
{
  psi->replaceGlobalValue ( p[k],
			    0,
			    psiView[k]
			  );
}
// Create the interface between NOX and the application
// This object is derived from NOX::Epetra::Interface
glsystem = Teuchos::rcp ( new GlSystem ( glProblem, eComm, reverse, psi ) );
}
else
glsystem = Teuchos::rcp ( new GlSystem ( glProblem, eComm, reverse ) );
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// setting the LOCA parameters
// ---------------------------------------------------------------------------
// Create LOCA sublist
Teuchos::RCP<Teuchos::ParameterList> paramList =
			     Teuchos::rcp(new Teuchos::ParameterList);
Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
  stepperList.set("Continuation Method", "Arc Length");// Default
  //stepperList.set("Continuation Method", "Natural");
  stepperList.set("Continuation Parameter", "H0");  // Must set
  stepperList.set("Initial Value", problemParameters.get<double>("H0"));     // Must set
  stepperList.set("Max Value",  3.0);                // Must set
  stepperList.set("Min Value", -3.0);                // Must set
  stepperList.set("Max Steps", 10000);               // Should set
  stepperList.set("Max Nonlinear Iterations", 20);   // Should set
  // ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// Setting the bifurcation parameters
// ---------------------------------------------------------------------------
// Create bifurcation sublist
Teuchos::ParameterList& bifurcationList =
				  locaParamsList.sublist("Bifurcation");
bifurcationList.set("Type", "None");                 // Default
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// Create predictor sublist
// ---------------------------------------------------------------------------
Teuchos::ParameterList& predictorList =
locaParamsList.sublist("Predictor");
//predictorList.set("Method", "Secant");               // Default
//predictorList.set("Method", "Constant");
predictorList.set("Method", "Tangent");
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// Create step size sublist
// ---------------------------------------------------------------------------
Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
stepSizeList.set("Method", "Adaptive");
double initialStepSize = 1.0e-4;
if (reverse)
stepSizeList.set("Initial Step Size", -initialStepSize);
else
stepSizeList.set("Initial Step Size",  initialStepSize);
stepSizeList.set("Min Step Size", 1.0e-4);
  stepSizeList.set("Max Step Size", 1.0e-3);
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Setting the nonlinear solver parameters
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

  Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
  if (verbose) {
      nlPrintParams.set( "Output Information",
                         NOX::Utils::Details +
                         NOX::Utils::OuterIteration +
                         NOX::Utils::InnerIteration +
                         NOX::Utils::Warning +
                         NOX::Utils::StepperIteration +
                         NOX::Utils::StepperDetails +
                         NOX::Utils::StepperParameters );  // Should set
  } else {
      nlPrintParams.set( "Output Information", 0 );  // No output of Trilinos
  }

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Sublist for linear solver for the Newton method
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  // lsParams.set("Aztec Solver", "BiCGStab");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Output Frequency", 32);
  lsParams.set("Output Solver Details", true);
  lsParams.set("Max Iterations", 2000);
  lsParams.set("Tolerance", 1e-4);

  // lsParams.set("Use Preconditioner as Solver",true);
  lsParams.set("Max Age Of Prec", 1);

  lsParams.set("Preconditioner","New Ifpack");
  lsParams.set("Ifpack Preconditioner","Amesos");
  Teuchos::ParameterList& IFPACKparams =  lsParams.sublist("Ifpack");
  IFPACKparams.set("amesos: solver type","Amesos_Superlu");
  // IFPACKparams.set("fact: level-of-fill",1);
  IFPACKparams.set("fact: ilut level-of-fill",1.0);
  IFPACKparams.set("fact: drop tolerance", 1e-1);
  IFPACKparams.set("schwarz: combine mode", "Zero");
  IFPACKparams.set("schwarz: compute condest", true);
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Let's force all status tests to do a full check
  nlParams.sublist("Solver Options").set("Status Test Check Type", "Complete");
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Create the necessary objects
  // ---------------------------------------------------------------------------
  // Create Epetra factory
  Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
                                        Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  Teuchos::RCP<LOCA::GlobalData> globalData =
                               LOCA::createGlobalData(paramList, epetraFactory);

  // set the directory to which all output gets written
  glsystem->setOutputDir( outputdir );

  // get the initial solution
  Teuchos::RCP<Epetra_Vector> soln = glsystem->getSolution();
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
#ifdef HAVE_LOCA_ANASAZI
  if (computeEigenvalues) {
    // Create Anasazi Eigensolver sublist
    stepperList.set("Compute Eigenvalues",true);
    Teuchos::ParameterList& aList = stepperList.sublist("Eigensolver");
    aList.set("Method", "Anasazi");
    if (!verbose)
      aList.set("Verbosity", Anasazi::Errors);

//    aList.set("Operator","Cayley");
//    double sigma = 10.0;
 //   aList.set("Cayley Pole",sigma);
//    double mu = 11.0;
//    aList.set("Cayley Zero",mu);

    aList.set("Operator","Jacobian Inverse");

    aList.set("Num Eigenvalues", 20);
    aList.set("Sorting Order", "LM"); // largest magnitude
    aList.set("Num Blocks", 30 ); // = max # of Arnoldi steps
    aList.set("Maximum Restarts", 500);

    aList.set("Save Eigen Data Method","User-Defined");
    aList.set("User-Defined Save Eigen Data Name", "MySave");

//     Teuchos::RCP<LOCA::Parameter::SublistParser> parser =
//         Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
//    Teuchos::RCP<Teuchos::ParameterList> aListPtr = Teuchos::rcp(&aList,false);
//

    Teuchos::RCP<Teuchos::ParameterList> stepperListPtr = 
        Teuchos::RCP<Teuchos::ParameterList>( &stepperList );
        //Teuchos::RCP<Teuchos::ParameterList>( new Teuchos::ParameterList(stepperList) );
//cout << "IN MAIN " << endl;
//cout << *stepperListPtr;
//cout << "END MAIN " << endl;
    std::string fileName = outputdir + "/eigenvalues.dat";
    Teuchos::RCP<EigenSaver> yourGreatSaver =
                           Teuchos::RCP<EigenSaver>( new EigenSaver(stepperListPtr,globalData,fileName,glsystem) );

    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> myGreatSaver =
                                                                 yourGreatSaver;
    aList.set("MySave",myGreatSaver);
  }
#else
    stepperList.set("Compute Eigenvalues",false);
#endif
  
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Create all possible Epetra_Operators.
  Teuchos::RCP<Epetra_RowMatrix> Analytic = glsystem->getJacobian();

  // Create the linear system
  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = glsystem;
  Teuchos::RCP<NOX ::Epetra::Interface::Jacobian> iJac = glsystem;

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys
    = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams,
                                                        lsParams,
                                                        glsystem,
                                                        iJac,
                                                        Analytic,
                                                        *soln )
                  );
  // ---------------------------------------------------------------------------



  // ---------------------------------------------------------------------------
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  LOCA::ParameterVector locaParams;
  locaParams.addParameter( "H0", problemParameters.get<double>("H0") );

  NOX::Epetra::Vector initialGuess(soln, NOX::Epetra::Vector::CreateView);

  Teuchos::RCP<LOCA::Epetra::Group> grp =
    Teuchos::rcp(new LOCA::Epetra::Group( globalData,
                                          nlPrintParams,
                                          iReq,
                                          initialGuess,
                                          linSys,
                                          locaParams )
                );

  grp->setParams( locaParams );
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Get the vector from the Problem
  Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
    Teuchos::rcp(new NOX::Epetra::Vector(soln,
                                         NOX::Epetra::Vector::CreateView));
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Set up the status tests
  // ---------------------------------------------------------------------------
  Teuchos::RCP<NOX::StatusTest::NormF> normF =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-13));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
  Teuchos::RCP<NOX::StatusTest::Generic> comboOR =
    Teuchos::rcp(new NOX::StatusTest::Combo( NOX::StatusTest::Combo::OR,
                                             normF,
                                             maxIters));
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Create the stepper
  LOCA::Stepper stepper(globalData, grp, comboOR, paramList);
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Perform continuation run
  LOCA::Abstract::Iterator::IteratorStatus status;
  try {
     status = stepper.run();
  }
  catch ( char const* e ) {
     std::cerr << "Exception raised: " << e << std::endl;
  }
  // ---------------------------------------------------------------------------

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
