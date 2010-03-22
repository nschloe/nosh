#include "Ginla_IO_SaveEigenData.h"

#include <NOX_Abstract_MultiVector.H>

// =============================================================================
Ginla::IO::SaveEigenData::
SaveEigenData ( Teuchos::RCP<Teuchos::ParameterList>           & eigenParamList,
                const Teuchos::RCP<const Recti::Grid::General> & grid,
                const Teuchos::RCP<const Ginla::Komplex>       & komplex,
                const Teuchos::RCP<Ginla::IO::StatsWriter>     & statsWriter,
                const Teuchos::RCP<Ginla::IO::StateWriter>     & stateWriter
              ) :
        eigenParamList_ ( eigenParamList ),
        grid_ ( grid ),
        komplex_ ( komplex ),
        statsWriter_ ( statsWriter ),
        stateWriter_ ( stateWriter ),
        locaStepper_ ( Teuchos::null ),
        numComputeStableEigenvalues_ ( 3 ),
        maxEigenvaluesSave_ ( 20 )
{
}
// =============================================================================
Ginla::IO::SaveEigenData::~SaveEigenData()
{
}
// =============================================================================
void
Ginla::IO::SaveEigenData::
setLocaStepper ( const Teuchos::RCP<LOCA::Stepper> locaStepper )
{
    locaStepper_ = locaStepper;
}
// =============================================================================
void
Ginla::IO::SaveEigenData::
releaseLocaStepper()
{
    locaStepper_ = Teuchos::null;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
Ginla::IO::SaveEigenData::
save ( Teuchos::RCP<std::vector<double> >       & evals_r,
       Teuchos::RCP<std::vector<double> >       & evals_i,
       Teuchos::RCP<NOX::Abstract::MultiVector> & evecs_r,
       Teuchos::RCP<NOX::Abstract::MultiVector> & evecs_i
     )
{    
    unsigned int step = locaStepper_->getStepNumber();

    unsigned int numEigenValues = evals_r->size();

    // store the unstable eigenstate into files
    int numUnstableEigenvalues = 0;
    for ( unsigned int k = 0; k < numEigenValues; k++ )
    {
        if ( ( *evals_r ) [k] > 0.0 )
        {
            numUnstableEigenvalues++;

            // make sure that the imaginary part of the eigenvector
            // is in fact 0
            Teuchos::RCP<NOX::Abstract::Vector> imagPart =
                Teuchos::rcpFromRef ( ( *evecs_i ) [k] );
            TEUCHOS_ASSERT_EQUALITY( 0.0, imagPart->norm() );

            // transform the real part of the eigenvector into psi
            Teuchos::RCP<NOX::Abstract::Vector> realPart =
                Teuchos::rcpFromRef ( ( *evecs_r ) [k] );
            Teuchos::RCP<NOX::Epetra::Vector> realPartE =
                Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector> ( realPart, true );
            Teuchos::RCP<ComplexVector> psi =
                komplex_->real2complex( realPartE->getEpetraVector() );
                                                     
            stringstream eigenstateFileNameAppendix;
            eigenstateFileNameAppendix << "-eigenvalue" << numEigenValues;
            stateWriter_->write( psi,
                                 grid_,
                                 step,
                                 eigenstateFileNameAppendix.str() );
        }
    }

    // Create Teuchos::ParameterList containing the data to be put into the
    // stats file.
    Teuchos::RCP<Teuchos::ParameterList> eigenvaluesList;
    eigenvaluesList->set( "0Step", step );
    eigenvaluesList->set( "1#unstable ev", numUnstableEigenvalues );
    for ( unsigned int k = 0; k < maxEigenvaluesSave_; k++ )
    {
        std::stringstream label;
        label << "Re(lambda_" << k << ")";
        eigenvaluesList->set( label.str(), ( *evals_r ) [k] );
        
        label << "Im(lambda_" << k << ")";
        eigenvaluesList->set( label.str(), ( *evals_i ) [k] );
    }
    statsWriter_->setList( eigenvaluesList );
    statsWriter_->print();

//     eigenFileStream << step << "\t";
//     eigenFileStream << numUnstableEigenvalues << "\t";
// 
//     // Set the output format
//     // TODO Think about replacing this with NOX::Utils::Sci.
//     eigenFileStream.setf ( std::ios::scientific );
//     eigenFileStream.precision ( 15 );
// 
//     for ( unsigned int k = 0; k < min ( numEigenValues,maxEigenvaluesSave_ ); k++ )
//         eigenFileStream << "\t" << ( *evals_r ) [k] << "\t" << ( *evals_i ) [k];
// 
//     // print "NaN" as fill-ins if there are more columns than eigenvalues
//     if ( maxEigenvaluesSave_>numEigenValues )
//         for ( unsigned int k = 0; k < 2* ( maxEigenvaluesSave_-numEigenValues ); k++ )
//             eigenFileStream << "\tNaN                   ";
// 
//     eigenFileStream << std::endl;
//     eigenFileStream.close();

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
