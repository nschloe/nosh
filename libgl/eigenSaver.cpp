#include "eigenSaver.h"

#include <vector>

#include <NOX_Abstract_Group.H>
#include <NOX_Abstract_MultiVector.H>
#include <NOX_Utils.H>
#include <LOCA_GlobalData.H>

#include <EpetraExt_Utils.h>

// =============================================================================
EigenSaver::EigenSaver(const Teuchos::RCP<Teuchos::ParameterList> eigenParamList,
		               const Teuchos::RCP<LOCA::GlobalData>& globalData,
		               const std::string outputDir,
		               const std::string eigenvaluesFileName,
		               const std::string contFileBaseName,
		               const std::string eigenstateFileNameAppendix,
		               const Teuchos::RCP<GlSystem> glSys) :
eigenParamList_(eigenParamList),
outputDir_(outputDir),
eigenvaluesFilePath_(outputDir + "/" + eigenvaluesFileName),
contFileBaseName_(contFileBaseName),
eigenstateFileNameAppendix_(eigenstateFileNameAppendix),
globalData_(globalData),
glSys_(glSys),
locaStepper_(0),
numComputeStableEigenvalues_(3)
{
}
// =============================================================================
EigenSaver::~EigenSaver()
{
}
// =============================================================================
void
EigenSaver::setLocaStepper( const Teuchos::RCP<LOCA::Stepper> locaStepper )
{
	locaStepper_ = locaStepper;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
EigenSaver::save(Teuchos::RCP<std::vector<double> > &evals_r,
		         Teuchos::RCP<std::vector<double> > &evals_i,
		         Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_r,
		         Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_i
		         )
{
	// Keep track of how often this method is called.
	// This is actually somewhat ugly as it assumes that this number coincides
	// with the number of steps in the continuation.
	// Generally, though, this will probably be true.
	static int step = 0;
	step++;

	int numEigenValues = evals_r->size();

	std::ofstream eigenFileStream;

	if (step == 1) {
		eigenFileStream.open(eigenvaluesFilePath_.c_str(), ios::trunc);
		eigenFileStream << "# Step" << "\t#unstable ev";
		eigenFileStream << "\tRe(lambda_0)" << "\t\tIm(lambda_0)";
		for (int k = 1; k < numEigenValues; k++) {
			eigenFileStream << "\t\tRe(lambda_" << k << ")" << "\t\tIm(lambda_"
					<< k << ")";
		}
		eigenFileStream << std::endl;
	} else {
		// just append the contents to the file
		eigenFileStream.open(eigenvaluesFilePath_.c_str(), ios::app);
	}

	int numUnstableEigenvalues = 0;
	for (int k = 0; k < numEigenValues; k++) {
		if ((*evals_r)[k] > 0.0) {
			numUnstableEigenvalues++;
			std::string eigenstateFilePath = outputDir_ + "/"
					+ contFileBaseName_ + EpetraExt::toString(step) + "-"
					+ eigenstateFileNameAppendix_ + EpetraExt::toString(
					numUnstableEigenvalues) + ".vtk";

			Teuchos::RCP<NOX::Abstract::Vector> abVec = Teuchos::rcpFromRef(
					(*evecs_r)[k]);
			Teuchos::RCP<NOX::Epetra::Vector> myVec =
					Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector>(abVec, true);

			glSys_->writeAbstractStateToFile(myVec->getEpetraVector(),eigenstateFilePath);
		}
	}

	eigenFileStream << step << "\t";
	eigenFileStream << numUnstableEigenvalues << "\t";

	// Set the output format
	// TODO: Think about replacing this with NOX::Utils::Sci.
	eigenFileStream.setf(std::ios::scientific);
	eigenFileStream.precision(15);

	for (int k = 0; k < numEigenValues; k++)
		eigenFileStream << "\t" << (*evals_r)[k] << "\t" << (*evals_i)[k];

	eigenFileStream << std::endl;
	eigenFileStream.close();

	// Adapt the computation for the next step.
	// Make sure that approximately \c numComputeStableEigenvalues_ stable eigenvalues
	// will be computed in the next step.
	eigenParamList_->set("Num Eigenvalues", numUnstableEigenvalues + numComputeStableEigenvalues_ );

	// reset the eigensolver to take notice of the new values
	locaStepper_->eigensolverReset();

	return NOX::Abstract::Group::Ok;
}
// =============================================================================
