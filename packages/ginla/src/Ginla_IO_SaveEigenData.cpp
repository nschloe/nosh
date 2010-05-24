/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2009--2010  Nico Schl\"omer

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

#include "Ginla_IO_SaveEigenData.h"

#include "Ginla_Helpers.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_StateTranslator.h"

#include <NOX_Abstract_MultiVector.H>

#include <AnasaziSortManager.hpp>

// =============================================================================
Ginla::IO::SaveEigenData::
SaveEigenData ( Teuchos::RCP<Teuchos::ParameterList>             & eigenParamList,
                const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator,
                const Teuchos::RCP<const Ginla::IO::StateWriter> & stateWriter,
                const Teuchos::RCP<Ginla::IO::StatsWriter>       & statsWriter
              ) :
        eigenParamList_ ( eigenParamList ),
        stateTranslator_ ( stateTranslator ),
        stateWriter_ ( stateWriter ),
        statsWriter_ ( statsWriter ),
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
    // Can't fetch step index now, so rely on the function's
    // being called exactly once per step.
    static unsigned int step = -1;
    if ( !locaStepper_.is_null() )
        step = locaStepper_->getStepNumber();
    else
        step++;

    unsigned int numEigenValues = evals_r->size();

    // store the unstable eigenstate into files
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
            eigenstateFileNameAppendix << "-eigenvalue" << k;
            Teuchos::RCP<Ginla::State> eigenstate = stateTranslator_->createState( realPartE->getEpetraVector() );
            stateWriter_->write( eigenstate,
                                 step,
                                 eigenstateFileNameAppendix.str() );
                                        
            // The matrix and the eigenvalue is supposedly purely real,
            // so the eigenvector's real and imaginary parts are eigenvectors
            // in their own right. Check here for the imaginary part,
            // and print it, too, if it's nonzero.
            Teuchos::RCP<NOX::Abstract::Vector> imagPart =
                Teuchos::rcpFromRef ( ( *evecs_i ) [k] );
            if ( imagPart->norm() > 1.0e-15 )
            {
                Teuchos::RCP<NOX::Epetra::Vector> imagPartE =
                    Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector> ( imagPart, true );
                eigenstateFileNameAppendix << "-im";
                Teuchos::RCP<Ginla::State> eigenstate = stateTranslator_->createState( imagPartE->getEpetraVector() );
                stateWriter_->write( eigenstate,
                                     step,
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

    if ( !locaStepper_.is_null() )
    {
        // Adapt the computation for the next step.
        // Make sure that approximately \c numComputeStableEigenvalues_ stable eigenvalues
        // will be computed in the next step.
        // TODO Remove +1. This was introduced as a mere work-around the exception
        // /home/nico/software/trilinos/dev/source/piro-observer-improvements/packages/anasazi/src/AnasaziBlockKrylovSchur.hpp:1247:
        // 
        // Throw number = 1
        // 
        // Throw test that evaluated to true: ritzIndex_[numRitzVecs_-1]==1
        // 
        // Anasazi::BlockKrylovSchur::computeRitzVectors(): the number of required Ritz vectors splits a complex conjugate pair.
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
    
        // preserve the sort manager
        Teuchos::RCP<Anasazi::SortManager<double> > d =
            eigenParamList_->get<Teuchos::RCP<Anasazi::SortManager<double> > >( "Sort Manager" );
        
        // TODO For some reason, the  call to eigensolverReset destroys teh "Sort Manager" entry.
        //      No idea why. This is a potentially serious bug in Trilinos.
        // reset the eigensolver to take notice of the new values
        locaStepper_->eigensolverReset ( eigenParamList_ );
        
        eigenParamList_->set( "Sort Manager", d );
    }

    return NOX::Abstract::Group::Ok;
}
// =============================================================================
