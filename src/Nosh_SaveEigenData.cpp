// @HEADER
//
//    Helper class for computation and I/O of eigenvalues and -vectors.
//    Copyright (C) 2009--2012  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
// =============================================================================
#include "Nosh_SaveEigenData.hpp"

#include "Nosh_ModelEvaluator.hpp"
#include "Nosh_CsvWriter.hpp"

#include <NOX_Abstract_MultiVector.H>
#include <AnasaziSortManager.hpp>

namespace Nosh {
// =============================================================================
SaveEigenData::
SaveEigenData(Teuchos::ParameterList &eigenParamList,
              const Teuchos::RCP<const Nosh::ModelEvaluator> &modelEval,
              const Teuchos::RCP<Nosh::CsvWriter> &csvWriter
              ) :
  eigenParamListPtr_( Teuchos::rcpFromRef<Teuchos::ParameterList>(
                        eigenParamList) ),
  modelEval_( modelEval ),
  csvWriter_( csvWriter ),
  locaStepper_( Teuchos::null ),
  numComputeStableEigenvalues_( 6 )
{
}
// =============================================================================
SaveEigenData::
~SaveEigenData()
{
}
// =============================================================================
void
SaveEigenData::
setLocaStepper( const Teuchos::RCP<LOCA::Stepper> locaStepper )
{
  locaStepper_ = locaStepper;
}
// =============================================================================
void
SaveEigenData::
releaseLocaStepper()
{
  locaStepper_ = Teuchos::null;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
SaveEigenData::
save( Teuchos::RCP<std::vector<double> > &evals_r,
      Teuchos::RCP<std::vector<double> > &evals_i,
      Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_r,
      Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_i
      )
{
  // Can't fetch step index now, so rely on the function
  // being called exactly once per step.
  // Step number updated at the end of the function.
  static unsigned int step = 0;
  unsigned int numEigenValues = evals_r->size();

  // Consider eigenvalue above tol to be unstable, and between -tol and tol to be nullstable.
  // This is really loose here, and this makes sure that we get most actual nullvalues.
  // Sometimes they are really hard to approach using regular continuation, but if tol
  // flags an approximate nullvector, turning point continuation may help nailing it down.
  const double tol = 1.0e-5;

  // Store all eigenstates in files.
  unsigned int numStableEigenvalues = 0;
  unsigned int numUnstableEigenvalues = 0;
  unsigned int numNullvalues = 0;
  for ( unsigned int k = 0; k < numEigenValues; k++ )
  {
    double eigenvalue = (*evals_r) [k];
    stringstream eigenstateFileNameAppendix;
    if ( eigenvalue  < -tol )
      // don't consider stable values
      eigenstateFileNameAppendix << "-seigenstate" << numStableEigenvalues++;
    else if ( eigenvalue  > tol )
      eigenstateFileNameAppendix << "-ueigenstate" << numUnstableEigenvalues++;
    else
      eigenstateFileNameAppendix << "-nullstate" << numNullvalues++;

    // transform the real part of the eigenvector into psi
    Teuchos::RCP<NOX::Abstract::Vector> realPart =
      Teuchos::rcpFromRef( (*evecs_r) [k] );
    Teuchos::RCP<NOX::Epetra::Vector> realPartE =
      Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector> ( realPart, true );

    // Don't store eigendata in a file right now.
    // The reason for this is that that the data gets written to the same
    // file as the solution because the mesh (sitting in the model
    // evaluator) determines the file name.
    // To fix this, there are several options:
    //   1) Store the eigen data in the same file, but into different
    //      fields. To this end, besides "psi_", several other fields would
    //      have to be declared and filled, e.g., "eig0_", "eig1_", ...
    //   2) Put the eigenstates into different files. This would probably
    //      require maintaining different mesh instances.
    //   3) A mix of the two first options: Keep the solution states
    //      separate and put all the eigenstates into one file.
    //
    //Teuchos::RCP<Nosh::State> eigenstate =
    //    modelEval_->createSavable( realPartE->getEpetraVector() );
    //eigenstate->save( step );  //  eigenstateFileNameAppendix.str();

    // The matrix and the eigenvalue is supposedly purely real,
    // so the eigenvector's real and imaginary parts are eigenvectors
    // in their own right. Check here for the imaginary part,
    // and print it, too, if it's nonzero.
    Teuchos::RCP<NOX::Abstract::Vector> imagPart =
      Teuchos::rcpFromRef( (*evecs_i) [k] );
    TEUCHOS_ASSERT_INEQUALITY( imagPart->norm(), <, 1.0e-15 );
//        if ( imagPart->norm() > 1.0e-15 )
//        {
//            Teuchos::RCP<NOX::Epetra::Vector> imagPartE =
//                Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector> ( imagPart, true );
//            eigenstateFileNameAppendix << "-im";
//            Teuchos::RCP<Nosh::State> eigenstate =
//                modelEval_->createSavable( imagPartE->getEpetraVector() );
//            eigenstate->save( step ); //eigenstateFileNameAppendix.str()
//        }
  }

  // Create Teuchos::ParameterList containing the data to be put into the
  // stats file.
  Teuchos::ParameterList eigenvaluesList;
  eigenvaluesList.set( "##step", step );
  eigenvaluesList.set( "#0unstable", numUnstableEigenvalues );
  eigenvaluesList.set( "#1null", numNullvalues );
  eigenvaluesList.set( "#2stable", numStableEigenvalues );
  for ( unsigned int k = 0; k < numEigenValues; k++ )
  {
    std::stringstream label;
    label << setw( 2 ) << setfill( '0' ) << k << "-0Re()";
    eigenvaluesList.set( label.str(), (*evals_r) [k] );

    // make sure that the imaginary part is indeed 0
    TEUCHOS_ASSERT_INEQUALITY( fabs( (*evals_i) [k] ), <, 1.0e-15 );
  }

  // Write out the data.
  if (step == 0)
    csvWriter_->writeHeader(eigenvaluesList);
  csvWriter_->writeRow(eigenvaluesList);

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

  // Make sure that the shift SIGMA (if using Shift-Invert) sits THRESHOLD above
  // the rightmost eigenvalue.
  std::string &op = eigenParamListPtr_->get<string> ( "Operator" );
  if ( !locaStepper_.is_null() && op.compare( "Shift-Invert" ) == 0 )
  {
    double maxEigenval = *std::max_element( evals_r->begin(), evals_r->end() );
    double threshold = 0.5;
    eigenParamListPtr_->set( "Shift", maxEigenval + threshold );

    // Preserve the sort manager.
    // TODO For some reason, the call to eigensolverReset destroys the "Sort Manager" entry.
    //      No idea why. This is potentially a bug in Trilinos.
    Teuchos::RCP<Anasazi::SortManager<double> > d =
      eigenParamListPtr_->get<Teuchos::RCP<Anasazi::SortManager<double> > >(
        "Sort Manager" );
    // reset the eigensolver to take notice of the new values
    locaStepper_->eigensolverReset( eigenParamListPtr_ );
    eigenParamListPtr_->set( "Sort Manager", d );
  }

  // update the step index for the next run
  if ( !locaStepper_.is_null() )
    step = locaStepper_->getStepNumber();
  else
    step++;

  return NOX::Abstract::Group::Ok;
}
// =============================================================================
} // namespace Nosh