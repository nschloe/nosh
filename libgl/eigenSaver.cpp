#include "eigenSaver.h"

#include <vector>

#include <NOX_Abstract_Group.H>
#include <NOX_Abstract_MultiVector.H>

#include <EpetraExt_Utils.h>

// =============================================================================
EigenSaver::EigenSaver( const std::string fileName ) :
  fileName_(fileName)
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
  for (int k=0; k<numEigenValues; k++ )
      if ( (*evals_r)[k]>0.0 )
          numUnstableEigenvalues++;

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

  return NOX::Abstract::Group::Ok;
}
// =============================================================================