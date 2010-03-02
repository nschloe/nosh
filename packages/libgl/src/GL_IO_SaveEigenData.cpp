#include "GL_IO_SaveEigenData.h"

#include <vector>

#include <NOX_Abstract_Group.H>
#include <NOX_Abstract_MultiVector.H>
#include <NOX_Utils.H>
#include <LOCA_GlobalData.H>

// =============================================================================
GL::IO::SaveEigenData::SaveEigenData ( const Teuchos::RCP<Teuchos::ParameterList> eigenParamList,
                         const std::string  outputDir,
                         const std::string  eigenvaluesFileName,
                         const std::string  contFileBaseName,
                         const std::string  eigenstateFileNameAppendix,
                         const Teuchos::RCP<AbstractStateWriter> stateWriter ) :
        eigenParamList_ ( eigenParamList ),
        outputDir_ ( outputDir ),
        eigenvaluesFilePath_ ( outputDir + "/" + eigenvaluesFileName ),
        contFileBaseName_ ( contFileBaseName ),
        eigenstateFileNameAppendix_ ( eigenstateFileNameAppendix ),
        stateWriter_ ( stateWriter ),
        locaStepper_ ( 0 ),
        numComputeStableEigenvalues_ ( 3 ),
        maxEigenvaluesSave_ ( 20 )
{
}
// =============================================================================
GL::IO::SaveEigenData::~SaveEigenData()
{
}
// =============================================================================
void
GL::IO::SaveEigenData::setLocaStepper ( const Teuchos::RCP<LOCA::Stepper> locaStepper )
{
    locaStepper_ = locaStepper;
}
// =============================================================================
void
GL::IO::SaveEigenData::releaseLocaStepper()
{
    locaStepper_ = Teuchos::null;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
GL::IO::SaveEigenData::save ( Teuchos::RCP<std::vector<double> > &evals_r,
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

    unsigned int numEigenValues = evals_r->size();

    std::ofstream eigenFileStream;

    if ( step == 1 )
    {
        eigenFileStream.open ( eigenvaluesFilePath_.c_str(), ios::trunc );
        eigenFileStream << "# Step" << "\t#unstable ev";
        for ( unsigned int k = 0; k < maxEigenvaluesSave_; k++ )
        {
            if ( k==0 )
                eigenFileStream << "\t";
            else
                eigenFileStream << "\t\t";
            eigenFileStream << "Re(lambda_" << k << ")" << "\t\tIm(lambda_"
            << k << ")";
        }
        eigenFileStream << std::endl;
    }
    else
    {
        // just append the contents to the file
        eigenFileStream.open ( eigenvaluesFilePath_.c_str(), ios::app );
    }

    int numUnstableEigenvalues = 0;
    for ( unsigned int k = 0; k < numEigenValues; k++ )
    {
        if ( ( *evals_r ) [k] > 0.0 )
        {
            numUnstableEigenvalues++;
            stringstream eigenstateString;
            eigenstateString
            << outputDir_  << "/"
            << contFileBaseName_ << step << "-"
            << eigenstateFileNameAppendix_ 
            << numUnstableEigenvalues << ".vtk";

            Teuchos::RCP<NOX::Abstract::Vector> abVec =
                Teuchos::rcpFromRef ( ( *evecs_r ) [k] );
            Teuchos::RCP<NOX::Epetra::Vector> myVec =
                Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector> ( abVec, true );

            stateWriter_->writeAbstractStateToFile ( myVec->getEpetraVector(),
                                                     eigenstateString.str() );
        }
    }

    eigenFileStream << step << "\t";
    eigenFileStream << numUnstableEigenvalues << "\t";

    // Set the output format
    // TODO: Think about replacing this with NOX::Utils::Sci.
    eigenFileStream.setf ( std::ios::scientific );
    eigenFileStream.precision ( 15 );

    for ( unsigned int k = 0; k < min ( numEigenValues,maxEigenvaluesSave_ ); k++ )
        eigenFileStream << "\t" << ( *evals_r ) [k] << "\t" << ( *evals_i ) [k];

    // print "NaN" as fill-ins if there are more columns than eigenvalues
    if ( maxEigenvaluesSave_>numEigenValues )
        for ( unsigned int k = 0; k < 2* ( maxEigenvaluesSave_-numEigenValues ); k++ )
            eigenFileStream << "\tNaN                   ";

    eigenFileStream << std::endl;
    eigenFileStream.close();

    // Adapt the computation for the next step.
    // Make sure that approximately \c numComputeStableEigenvalues_ stable eigenvalues
    // will be computed in the next step.
    int nextNumEigenvalues = numUnstableEigenvalues + numComputeStableEigenvalues_;
    eigenParamList_->set ( "Num Eigenvalues", nextNumEigenvalues );

    // Make sure that the shift SIGMA (if using Shift-Invert) sits THRESHOLD above
    // the rightmost eigenvalue.
    if ( eigenParamList_->get<string> ( "Operator" ).compare ( "Shift-Invert" ) ==0 )
    {
        double maxEigenval = *std::max_element ( evals_r->begin(), evals_r->end() );
        double threshold = 0.5;
        eigenParamList_->set ( "Shift", maxEigenval + threshold );
    }

    // reset the eigensolver to take notice of the new values
    locaStepper_->eigensolverReset ( eigenParamList_ );

    return NOX::Abstract::Group::Ok;
}
// =============================================================================
