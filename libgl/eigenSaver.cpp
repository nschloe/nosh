#include "eigenSaver.h"

#include <vector>

#include <NOX_Abstract_Group.H>
#include <NOX_Abstract_MultiVector.H>
#include <NOX_Utils.H>
#include <LOCA_GlobalData.H>

#include <EpetraExt_Utils.h>

// =============================================================================
EigenSaver::EigenSaver( const Teuchos::ParameterList& eigenParamsList,
                        const Teuchos::RCP<LOCA::GlobalData>& globalData,
		        const std::string fileName,
                        const Teuchos::RCP<GlSystem> glSys ) :
  eigenParamsList_(eigenParamsList),
  fileName_(fileName),
  globalData_(globalData),
  glSys_(glSys)
  //parsedParams_(Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData))),
{
};
// =============================================================================
EigenSaver::~EigenSaver()
{
};
// =============================================================================
NOX::Abstract::Group::ReturnType
EigenSaver::save ( Teuchos::RCP<std::vector<double> >       &evals_r,
                   Teuchos::RCP<std::vector<double> >       &evals_i,
                   Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_r,
                   Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_i  )
{

// std::cout << "Hi before" << std::endl;
// std::cout << globalData_->locaUtils->StepperParameters;
// globalData_->locaUtils->print(std::cout);
// std::cout << "Hi after" << std::endl;
//globalData_->locaUtils();//->StepperParameters();

  // Keep track of how often this method is called.
  // This is actually somewhat ugly as it assumes that this number coincides
  // with the number of steps in the continuation.
  // Generally, though, this will probably be true.
  static int step = 0;
  step++;

  int numEigenValues = evals_r->size();

  std::ofstream eigenFileStream;

  if ( abs(step)==1 ) {
      eigenFileStream.open (fileName_.c_str(),ios::trunc);
      eigenFileStream << "# Step"
                      << "\t#unstable ev";
      eigenFileStream <<   "\tRe(lambda_0)"
                      << "\t\tIm(lambda_0)";
      for (int k=1; k<numEigenValues; k++ ) {
          eigenFileStream << "\t\tRe(lambda_" << k << ")"
                          << "\t\tIm(lambda_" << k << ")";
      }
      eigenFileStream << std::endl;
  } else {
      // just append to the the contents to the file
      eigenFileStream.open (fileName_.c_str(),ios::app);
  }

  int numUnstableEigenvalues = 0;
  for (int k=0; k<numEigenValues; k++ ) {
      if ( (*evals_r)[k]>0.0 ) {
          numUnstableEigenvalues++;
          std::string eigenStateFileName = "step-"
                                         + EpetraExt::toString(step)
                                         + "eigenfunction-"
                                         + EpetraExt::toString(numUnstableEigenvalues)
                                         + ".vtk";

          Teuchos::RCP<NOX::Abstract::Vector> abVec = Teuchos::rcpFromRef( (*evecs_r)[k] );
          Teuchos::RCP<NOX::Epetra::Vector> myVec = Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector>( abVec,true );

          Teuchos::RCP<Teuchos::ParameterList> tmpList;
          glSys_->printState( myVec->getEpetraVector(),
                              eigenStateFileName,
                              tmpList );

  }
  }

  eigenFileStream << step << "\t";
  eigenFileStream << numUnstableEigenvalues << "\t";

  // Set the output format
  // TODO: Think about replacing this with NOX::Utils::Sci.
  eigenFileStream.setf( std::ios::scientific );
  eigenFileStream.precision(15);

  for (int k=0; k<numEigenValues; k++ ) {
      eigenFileStream << "\t" << (*evals_r)[k]
                      << "\t" << (*evals_i)[k];
  }

  eigenFileStream << std::endl;
  eigenFileStream.close();

//  Teuchos::RCP<Teuchos::ParameterList> eigenParamsList = 
//	  Teuchos::RCP<Teuchos::ParameterList> ( new Teuchos::ParameterList( topLevelParams_->sublist("LOCA").
//								      sublist("Stepper",true).
//								      sublist("Eigensolver",true) ) ;
  
//  Teuchos::ParameterList & eigenParamsList = topLevelParams_->sublist("LOCA").
//					         	     sublist("Stepper",true).
//							     sublist("Eigensolver",true);
  cout << eigenParamsList_ ;
  //eigenParamsList.set("Num Eigenvalues",step);
  
  return NOX::Abstract::Group::Ok;
}
// =============================================================================
// void
// EigenSaver::saveEigenstate ( const std::string                         fileName,
//                              const Teuchos::RCP<NOX::Abstract::Vector> &evec_r  )
// {
// //   conParam = 0.0;
// //   glSys->GlSystem::printSolution ( evec_r, conParam );
// 
//   // create complex vector
//   Teuchos::RCP<const Tpetra::Map<int> > ComplexMap = glSys_->getComplexMap();
//   Tpetra::MultiVector<double_complex,int>  psi(ComplexMap,1,true);
// 
//   glSys->real2complex ( evec_r, psi );
// 
//   // create parameter list to be written to the file
//   Teuchos::ParameterList tmpList;
// //   tmpList.get ( "edgelength", glSys_->getStaggeredGrid()->getEdgeLength() );
// //   tmpList.get ( "Nx",         Gl_.getStaggeredGrid()->getNx() );
// 
//   IoVirtual* fileIo = IoFactory::createFileIo ( outputDir_+"/"+fileName );
//   fileIo->write ( psi,
//                   tmpList,
//                   * ( Gl_.getStaggeredGrid() ) );
// 
// }
// // =============================================================================
