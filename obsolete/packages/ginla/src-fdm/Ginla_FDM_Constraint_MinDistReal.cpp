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

#include "Ginla_Constraint_MinDistReal.h"

#include "Ginla_LocaSystem_Default.h"

#include <NOX_Epetra_MultiVector.H>

// =============================================================================
// Default constructor
Ginla::Constraint::MinDistReal::
MinDistReal ( const Teuchos::RCP<const Epetra_Vector>  & x,
              const LOCA::ParameterVector              & paramsVector
            ):
  constraints_(1,1),
  isValidConstraints_(false),
  paramsVector_( paramsVector ),
  // Set the initial guess.
  // psiRef_ is adapted accordingly before the first continuation step
  // (see preProcessContinuationStep()).
  x_( Teuchos::rcp( new NOX::Epetra::Vector( Epetra_Vector(*x) ) ) ),
  xRef_( Teuchos::rcp( new NOX::Epetra::Vector( Epetra_Vector(x->Map()) ) ) ),
  xRefDot_( Teuchos::rcp( new NOX::Epetra::MultiVector( Epetra_Vector(x->Map(),1) ) ) )
{ 
  // Initialize constraint
  constraints_.putScalar(0.0);
  return;
}
// =============================================================================
// Destructor
Ginla::Constraint::MinDistReal::
~MinDistReal()
{
}
// =============================================================================
// TODO what does this function do? deep copy?
Ginla::Constraint::MinDistReal::
MinDistReal( const Ginla::Constraint::MinDistReal & source,
             NOX::CopyType                          type
           ) :
  constraints_(source.constraints_),
  isValidConstraints_(false),
  paramsVector_( source.paramsVector_ ),
  x_( Teuchos::rcp( new NOX::Abstract::Vector(*source.x_) ) ),
  xRef_( Teuchos::rcp( new NOX::Epetra::Vector(*source.xRef_) ) ),
  xRefDot_( Teuchos::rcp( new NOX::Epetra::MultiVector(*source.xRefDot_) ) )
{
  if (source.isValidConstraints_ && type == NOX::DeepCopy)
    isValidConstraints_ = true;
  return;
}
// =============================================================================
void
Ginla::Constraint::MinDistReal::
copy( const LOCA::MultiContinuation::ConstraintInterface & src )
{ 
  const MinDistReal & source = Teuchos::dyn_cast<const MinDistReal>( src );

  if (this != &source)
  {
    // deep copies:
    constraints_ = source.constraints_;
    isValidConstraints_ = source.isValidConstraints_;
    paramsVector_ = source.paramsVector_;
    *x_ = *source.x_;
    *xRef_ = *source.xRef_;
    *xRefDot_ = *source.xRefDot_;
  }

  return;
}
// =============================================================================
Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
Ginla::Constraint::MinDistReal::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Ginla::Constraint::MinDistReal(*this, type));
}
// =============================================================================
void
Ginla::Constraint::MinDistReal::
preProcessContinuationStep( LOCA::Abstract::Iterator::StepStatus stepStatus )
{
  // deep copy the current guess as a reference solution
  *xRef_ = *x_;

  isValidDx_ = false;
  isValidConstraints_ = false;
  return;
}
// =============================================================================
int
Ginla::Constraint::MinDistReal::
numConstraints() const
{
  return constraints_.numRows();
}
// =============================================================================
void
Ginla::Constraint::MinDistReal::
setX(const NOX::Abstract::Vector & y)
{ 
  *x_ = y;
  return;
}
// =============================================================================
void
Ginla::Constraint::MinDistReal::
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
Ginla::Constraint::MinDistReal::
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
Ginla::Constraint::MinDistReal:: 
computeConstraints()
{
  if (!isValidConstraints_)
  {
    computeDX(); // update xRefDot_

    
    std::cout  << "111" << std::endl;
//     NOX::Epetra::Vector x = Teuchos::dyn_cast<NOX::Epetra::Vector>((*xRefDot_)[0]);
    xRefDot_->print( std::cout );
    std::cout  << "1b" << std::endl;
    constraints_(0,0) =  (*xRefDot_)[0].innerProduct( x_ );
    std::cout << "222" << std::endl;
    isValidConstraints_ = true;
  }

  return NOX::Abstract::Group::Ok;
}
// =============================================================================
const NOX::Abstract::MultiVector*
Ginla::Constraint::MinDistReal::
getDX () const
{
    return xRefDot_.get();
}
// =============================================================================
NOX::Abstract::Group::ReturnType
Ginla::Constraint::MinDistReal::
computeDX()
{
  if (!isValidDx_)
  {
      // \Im( psi0^H psi ) = [-Im(psi0),Re(psi0)] [ Re(psi); Im(psi) ].
   
      // Get the underlying Epetra_Vectors to gain access to the elements.
      Epetra_Vector & xRefE    = xRef_->getEpetraVector();
      Teuchos::RCP<Epetra_Vector> xRefDotE = Teuchos::rcp( (xRefDot_->getEpetraMultiVector())(0) );
      
      TEUCHOS_ASSERT( xRefE.GlobalLength() % 2 == 0 );

      int N = xRefE.GlobalLength()/2;
      for ( int k=0; k<N; k++ )
      {
          (*xRefDotE)[2*k]   = -xRefE[2*k+1]; // -Im(psi0)
          (*xRefDotE)[2*k+1] =  xRefE[2*k];   //  Re(psi_0)
      }

      isValidDx_ = true;
  }

  return NOX::Abstract::Group::Ok;
}
// =============================================================================
// The first column of dgdp should be filled with the constraint
// residuals g if isValidG is false.
// If isValidG is true, then the dgdp contains g on input.
NOX::Abstract::Group::ReturnType
Ginla::Constraint::MinDistReal::
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
Ginla::Constraint::MinDistReal::
isConstraints() const
{
  return isValidConstraints_;
}
// =============================================================================
bool
Ginla::Constraint::MinDistReal::
isDX() const
{
  return isValidDx_;
}
// =============================================================================
const NOX::Abstract::MultiVector::DenseMatrix &
Ginla::Constraint::MinDistReal::
getConstraints() const
{
  return constraints_;
}
// =============================================================================
// Return true if solution component of constraint derivatives is zero.
bool
Ginla::Constraint::MinDistReal::
isDXZero() const
{
  return false;
}
// =============================================================================
