/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "Ginla_BranchExplorer_ProcessEigenData.h"

#include <LOCA_Stepper.H>
#include <NOX_Abstract_MultiVector.H>
#include <NOX_Epetra_Vector.H>

#include "Ginla_IO_EigenSaver_Abstract.h"
#include "Ginla_Helpers.h"
#include "Ginla_IO_StatsWriter.h"

// =============================================================================
Ginla::BranchExplorer::ProcessEigenData::
ProcessEigenData ( Teuchos::RCP<Teuchos::ParameterList>                      & eigenParamList,
                   const Teuchos::RCP<const Ginla::IO::EigenSaver::Abstract> & eigenSaver,
                   const Teuchos::RCP<Ginla::IO::StatsWriter>                & statsWriter
                 ) :
        eigenParamList_ ( eigenParamList ),
        eigenSaver_ ( eigenSaver ),
        statsWriter_ ( statsWriter ),
        locaStepper_ ( Teuchos::null ),
        crossingEigenvectors_( Teuchos::Array<Epetra_Vector>() ),
        previousNumUnstableEigenvalues_( 0 ),
        numComputeStableEigenvalues_ ( 3 ),
        maxEigenvaluesSave_ ( 20 )
{
}
// =============================================================================
Ginla::BranchExplorer::ProcessEigenData::
~ProcessEigenData()
{
}
// =============================================================================
void
Ginla::BranchExplorer::ProcessEigenData::
setLocaStepper ( const Teuchos::RCP<LOCA::Stepper> locaStepper )
{
    locaStepper_ = locaStepper;
}
// =============================================================================
void
Ginla::BranchExplorer::ProcessEigenData::
releaseLocaStepper()
{
    locaStepper_ = Teuchos::null;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
Ginla::BranchExplorer::ProcessEigenData::
save ( Teuchos::RCP<std::vector<double> >       & evals_r,
       Teuchos::RCP<std::vector<double> >       & evals_i,
       Teuchos::RCP<NOX::Abstract::MultiVector> & evecs_r,
       Teuchos::RCP<NOX::Abstract::MultiVector> & evecs_i
     )
{    
    unsigned int step = locaStepper_->getStepNumber();

    unsigned int numEigenValues = evals_r->size();

    // store the unstable eigenstates into files
    unsigned int numUnstableEigenvalues = 0;
    for ( unsigned int k = 0; k < numEigenValues; k++ )
    {
        if ( ( *evals_r ) [k] > 0.0 )
        {
            numUnstableEigenvalues++;

            // transform the real part of the eigenvector into psi
            Teuchos::RCP<NOX::Abstract::Vector> realPart =
                Teuchos::rcpFromRef ( ( *evecs_r ) [k] );
            Teuchos::RCP<NOX::Epetra::Vector> realPartE =
                Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector> ( realPart, true );
                                                     
            stringstream eigenstateFileNameAppendix;
            eigenstateFileNameAppendix << "eigenvalue" << k;
            eigenSaver_->printSolution( realPartE->getEpetraVector(),
                                        eigenstateFileNameAppendix.str() );
                                        
            // The matrix and the eigenvalue is supposedly purely real,
            // so the eigenvector's real an imaginary parts are eigenvectors
            // in their own right. Check here for the imaginary part,
            // and print it, too, if it's nonzero.
            Teuchos::RCP<NOX::Abstract::Vector> imagPart =
                Teuchos::rcpFromRef ( ( *evecs_i ) [k] );
            if ( imagPart->norm() > 1.0e-15 )
            {
                Teuchos::RCP<NOX::Epetra::Vector> imagPartE =
                    Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector> ( imagPart, true );
                eigenstateFileNameAppendix << "-im";
                eigenSaver_->printSolution( imagPartE->getEpetraVector(),
                                            eigenstateFileNameAppendix.str() );
            }
        }
    }

    // Create Teuchos::ParameterList containing the data to be put into the
    // stats file.
    Teuchos::ParameterList eigenvaluesList;
    eigenvaluesList.set( "#step", step );
    eigenvaluesList.set( "#unstable", numUnstableEigenvalues );
    for ( unsigned int k = 0; k < maxEigenvaluesSave_; k++ )
    {
        std::stringstream label;
        label << setw( Ginla::Helpers::numDigits(maxEigenvaluesSave_) ) << setfill( '0' ) << k << "-0Re()";
        if ( k<numEigenValues )
            eigenvaluesList.set( label.str(), ( *evals_r ) [k] );
        else
            eigenvaluesList.set( label.str(), "----------------------" );
        
        // empty the stringstream
        label.str(std::string());
        
        label << setw( Ginla::Helpers::numDigits(maxEigenvaluesSave_) ) << setfill( '0' ) <<  k<< "-1Im()";
        if ( k<numEigenValues )
            eigenvaluesList.set( label.str(), ( *evals_i ) [k] );
        else
            eigenvaluesList.set( label.str(), "----------------------" );
    }
    statsWriter_->setList( eigenvaluesList );
    statsWriter_->print();
    

    // Check if eigenvalues crossed the imaginary axis. If yes, store
    // the corresponding eigenstates.
    unsigned int numCrossingEigenvalues;
    if ( numUnstableEigenvalues > previousNumUnstableEigenvalues_ )
    {
        // Find the abs(d) positive eigenvalues closest to the imaginary axis
        // and store the corresponding eigenstates.
        // Rely on the fact that the eigenvalues are in ascending order.
        numCrossingEigenvalues = numUnstableEigenvalues - previousNumUnstableEigenvalues_;
        unsigned int k0;
        for ( unsigned int k = 0; k < numEigenValues; k++ )
        {
            if ( ( *evals_r ) [k] > 0.0 )
            {
              k0 = k;
              break;
            }
        }
        // k0,...,k0+d-1 are the eigenvalues that have crossed.
        std::vector<unsigned int> indices( numCrossingEigenvalues );
        for ( unsigned int k = 0; k < numCrossingEigenvalues; k++ )
          indices[k] = k0 + k;
        
    }
    else if ( numUnstableEigenvalues < previousNumUnstableEigenvalues_ )
    {
        // Find the abs(d) negative eigenvalues closest to the imaginary axis
        // and store the corresponding eigenstates.
        // Rely on the fact that the eigenvalues are in ascending order.
        // TODO http://trilinos.sandia.gov/packages/docs/dev/packages/nox/doc/html/classAnasazi_1_1LOCASort.html
        numCrossingEigenvalues = previousNumUnstableEigenvalues_ - numUnstableEigenvalues;
        unsigned int k0;
        for ( unsigned int k = numEigenValues; k > 0; --k )
        {
            if ( ( *evals_r ) [k] < 0.0 )
            {
              k0 = k;
              break;
            }
        }
        // k0,...,k0-d+1 are the eigenvalues that have crossed.
        std::vector<unsigned int> indices( numCrossingEigenvalues );
        for ( unsigned int k = 0; k < numCrossingEigenvalues; k++ )
          indices[k] = k0 - k;
    }
    // store the eigenvectors 'indices'
    // store the unstable eigenstates into files
    crossingEigenvectors_.clear();
    for ( unsigned int k = 0; k < numCrossingEigenvalues; k++ )
    {
        Teuchos::RCP<NOX::Abstract::Vector> realPart =
            Teuchos::rcpFromRef ( ( *evecs_r ) [k] );
        Teuchos::RCP<NOX::Epetra::Vector> realPartE =
            Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector> ( realPart, true );

        crossingEigenvectors_.append( realPartE->getEpetraVector() );

        // Toss away the imaginary part:
        // Matrix and eigenvalue are purely real, hence the imaginary part,
        // just as well as the real part of the eigenvector, are
        // eigenvectors of the matrix in their own rights. We know that the
        // eigenspace as dimension d, and hoping that the real parts
        // of the eigenvectors considered here are linearly independent,
        // the imaginary parts wouldn't add anything.
    }

    // reset for the next step
    previousNumUnstableEigenvalues_ = numUnstableEigenvalues;


    // Adapt the computation for the next step.
    // Make sure that approximately \c numComputeStableEigenvalues_ stable eigenvalues
    // will be computed in the next step.
    int nextNumEigenvalues = numUnstableEigenvalues + numComputeStableEigenvalues_;
    eigenParamList_->set ( "Num Eigenvalues", nextNumEigenvalues );

    // Make sure that the shift SIGMA (if using Shift-Invert) sits THRESHOLD above
    // the rightmost eigenvalue.
    if ( eigenParamList_->get<string> ( "Operator" ).compare ( "Shift-Invert" ) == 0 )
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
Teuchos::Array<Epetra_Vector>
Ginla::BranchExplorer::ProcessEigenData::
getCrossingEigenvectors() const
{
    return crossingEigenvectors_; 
}
// =============================================================================