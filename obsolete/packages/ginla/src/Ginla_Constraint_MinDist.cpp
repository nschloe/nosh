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

#include "Ginla_Constraint_MinDist.h"

#include "Ginla_LocaSystem_Default.h"

#include <NOX_Epetra_MultiVector.H>

// =============================================================================
// Default constructor
Ginla::Constraint::MinDist::
MinDist ( const Teuchos::RCP<const Komplex2::LinearProblem> & komplex,
          const Teuchos::RCP<ComplexVector>                 & psi,
          const LOCA::ParameterVector                       & paramsVector
        ):
  komplex_( komplex ),
  constraints_(1,1),
  isValidConstraints_(false),
  // TODO: note that paramsVector must be a LOCA::ParameterVector containing all
  // parameters, and it should include chi
  paramsVector_( paramsVector ),
  // Set the initial guess.
  // psiRef_ is adapted accordingly before the first continuation step
  // (see preProcessContinuationStep()).
  psi_( psi ),
  psiRef_( Teuchos::rcp( new ComplexVector( psi->getMap() ) ) )
{ 
  // Initialize constraint
  constraints_.putScalar(0.0);
  return;
}
// =============================================================================
// Destructor
Ginla::Constraint::MinDist::
~MinDist()
{
}
// =============================================================================
// TODO what does this function do? deep copy?
Ginla::Constraint::MinDist::
MinDist( const Ginla::Constraint::MinDist & source,
         NOX::CopyType                      type
       ) :
  komplex_( Teuchos::rcp( new Komplex2::LinearProblem(*source.komplex_) ) ),
  constraints_(source.constraints_),
  isValidConstraints_(false),
  paramsVector_( source.paramsVector_ ),
  psi_( Teuchos::rcp( new ComplexVector(*source.psi_) ) ),
  psiRef_( Teuchos::rcp( new ComplexVector(*source.psiRef_) ) )
{
  if (source.isValidConstraints_ && type == NOX::DeepCopy)
    isValidConstraints_ = true;
  return;
}
// =============================================================================
void
Ginla::Constraint::MinDist::
copy( const LOCA::MultiContinuation::ConstraintInterface & src )
{ 
  const MinDist & source = Teuchos::dyn_cast<const MinDist>( src );

  if (this != &source)
  {
    // No deep copy needed as only
    //   komplex_->real2complex()
    // is used.
    komplex_ = source.komplex_;
    // deep copies:
    constraints_ = source.constraints_;
    isValidConstraints_ = source.isValidConstraints_;
    paramsVector_ = source.paramsVector_;
    *psi_ = *source.psi_;
    *psiRef_ = *source.psiRef_;
  }

  return;
}
// =============================================================================
Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
Ginla::Constraint::MinDist::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Ginla::Constraint::MinDist(*this, type));
}
// =============================================================================
void
Ginla::Constraint::MinDist::
preProcessContinuationStep( LOCA::Abstract::Iterator::StepStatus stepStatus )
{
  // deep copy the current guess as a reference solution
  *psiRef_ = *psi_;
  isValidConstraints_ = false;
  return;
}
// =============================================================================
int
Ginla::Constraint::MinDist::
numConstraints() const
{
  return constraints_.numRows();
}
// =============================================================================
void
Ginla::Constraint::MinDist::
setX(const NOX::Abstract::Vector & y)
{ 
  TEUCHOS_ASSERT( komplex_.is_valid_ptr() && !komplex_.is_null() );

  const Epetra_Vector & yE =
          Teuchos::dyn_cast<const NOX::Epetra::Vector>( y ).getEpetraVector();
  psi_ = komplex_->real2complex( yE );

  return;
}
// =============================================================================
void
Ginla::Constraint::MinDist::
setParam( int paramID, double val )
{
  paramsVector_[paramID] = val;
  
  // This would usually be set to FALSE, but in this case, the constraints
  // do not depend of the parameters. Hence, the constraints
  // remain valid even when changing the parameters.
  // isValidConstraints_ = false;
  
  return;
}
// =============================================================================
void
Ginla::Constraint::MinDist::
setParams( const vector<int>                             & paramIDs,
           const NOX::Abstract::MultiVector::DenseMatrix & vals )
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    paramsVector_[paramIDs[i]] = vals(i,0);

  // This would usually be set to FALSE, but in this case, the constraints
  // do not depend of the parameters. Hence, the constraints
  // remain valid even when changing the parameters.
  // isValidConstraints_ = false;
  
  return;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
Ginla::Constraint::MinDist::
computeConstraints()
{
  if (!isValidConstraints_)
  {
    // the dot() method is indeed smart enough to compute
    // conjugate the first complex vector psiRef_
    constraints_(0,0) = std::imag( psiRef_->dot(*psi_) );
    isValidConstraints_ = true;
  }

  return NOX::Abstract::Group::Ok;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
Ginla::Constraint::MinDist::
computeDX()
{
  return NOX::Abstract::Group::Ok;
}
// =============================================================================
// The first column of dgdp should be filled with the constraint
// residuals g if isValidG is false.
// If isValidG is true, then the dgdp contains g on input.
NOX::Abstract::Group::ReturnType
Ginla::Constraint::MinDist::
computeDP( const vector<int> & paramIDs,
           NOX::Abstract::MultiVector::DenseMatrix & dgdp, 
           bool isValidG
         )
{
  if (!isValidG)
      dgdp(0,0) = constraints_(0,0);

  for (unsigned int i=0; i<paramIDs.size(); i++)
      dgdp(0,i+1) = 0.0; // no dependence on constraints in the phase condition

  return NOX::Abstract::Group::Ok;
}
// =============================================================================
bool
Ginla::Constraint::MinDist::
isConstraints() const
{
  return isValidConstraints_;
}
// =============================================================================
bool
Ginla::Constraint::MinDist::
isDX() const
{
  return true;
}
// =============================================================================
const NOX::Abstract::MultiVector::DenseMatrix &
Ginla::Constraint::MinDist::
getConstraints() const
{
  return constraints_;
}
// =============================================================================
// Return true if solution component of constraint derivatives is zero.
bool
Ginla::Constraint::MinDist::
isDXZero() const
{
  return false;
}
// =============================================================================
// Compute result_p = alpha * dg/dx * input_x.
NOX::Abstract::Group::ReturnType
Ginla::Constraint::MinDist::
multiplyDX ( double                                    alpha,
             const NOX::Abstract::MultiVector        & input_x,
             NOX::Abstract::MultiVector::DenseMatrix & result_p
           ) const
{ 
  TEUCHOS_ASSERT( komplex_.is_valid_ptr() && !komplex_.is_null() );
  
  TEUCHOS_ASSERT_EQUALITY( result_p.numCols(), input_x.numVectors() );
  
  for ( int k=0; k<input_x.numVectors(); k++ )
  {
      const Epetra_Vector & xE =
              Teuchos::dyn_cast<const NOX::Epetra::Vector>( input_x[0] ).getEpetraVector();
      Teuchos::RCP<ComplexVector> xPsi = komplex_->real2complex( xE );

      result_p(0,k) = alpha * std::imag( psiRef_->dot(*xPsi) );
  }
  
  return NOX::Abstract::Group::Ok;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
Ginla::Constraint::MinDist::
addDX ( Teuchos::ETransp transb,
        double alpha,
        const NOX::Abstract::MultiVector::DenseMatrix &b,
        double beta,
        NOX::Abstract::MultiVector &result_x
      ) const
{
  TEST_FOR_EXCEPTION( true,
                      std::logic_error,
                      "Not yet implemented." );

  return NOX::Abstract::Group::Failed;
}
// =============================================================================
