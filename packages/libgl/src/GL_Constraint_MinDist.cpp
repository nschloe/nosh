/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

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

#include "GL_Constraint_MinDist.h"

#include "GL_LocaSystem_Default.h"

#include <NOX_Epetra_MultiVector.H>

// =============================================================================
// Default constructor
GL::Constraint::MinDist::MinDist ( const Teuchos::RCP <GL::LocaSystem::Default> & glSystem,
                                   const NOX::Abstract::Vector & initialGuess,
                                   const LOCA::ParameterVector & paramsVector ):
  glSystem_( glSystem ),
  constraints_(1,1),
  isValidConstraints_(false),
  // TODO: note that paramsVector must be a LOCA::ParameterVector containing all
  // parameters, and it should include chi
  paramsVector_( paramsVector ),
  psi_(Teuchos::null),
  psiRef_(Teuchos::null)
{ 
  // Initialize constraint
  constraints_.putScalar(0.0);

  // We're trying to get an Epetra_Vector here, but apparently
  // the straightforward cast NOX::Abstract::Vector-->NOX::Epetra::Vector
  // doesn't work. That works for MultiVector's though, so do the dance.
  NOX::Epetra::MultiVector & initialGuessE =
          dynamic_cast<NOX::Epetra::MultiVector&>(*(initialGuess.createMultiVector(1)));
  
  // At the first step, the reference solution is
  // the initial guess. For other steps, the solution
  // at the previous continuation step will be taken 
  // (see postProcessContinuationStep())
  psiRef_ = glSystem_->getGlKomplex()->real2complex( *(initialGuessE.getEpetraMultiVector()(0)) );
  
  // Soluion equals initialGuess at the instantiation.
  // Essentially we needed to shape x as xRef
  psi_ = glSystem_->getGlKomplex()->real2complex( *(initialGuessE.getEpetraMultiVector()(0)) );
}
// =============================================================================
// Destructor
GL::Constraint::MinDist::~MinDist()
{
}
// =============================================================================
// TODO what does this function do? deep copy?
GL::Constraint::MinDist::MinDist(const GL::Constraint::MinDist & source,
                                       NOX::CopyType             type ) :
  glSystem_( Teuchos::rcp( new GL::LocaSystem::Default(*source.glSystem_) ) ),
  constraints_(source.constraints_),
  isValidConstraints_(false),
  // TODO should we put paramsVector in here?
  psi_( Teuchos::rcp( new ComplexVector(*source.psi_) ) ),
  psiRef_( Teuchos::rcp( new ComplexVector(*source.psiRef_) ) )
{
  if (source.isValidConstraints_ && type == NOX::DeepCopy)
    isValidConstraints_ = true;
}
// =============================================================================
void
GL::Constraint::MinDist::copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  
  TEST_FOR_EXCEPT( true );
  
//   const GL::Constraint::MinDist& source =
//       dynamic_cast<const GL::Constraint::MinDist&>(src);
// 
//   if (this != &source) {
//     *glSystem_ = *source.glSystem_;
//     constraints_ = source.constraints_;
//     isValidConstraints_ = source.isValidConstraints_;
//     *psi_ = *source.psi_;
//     *psiRef_ = *source.psiRef_;
//   }

  return;
}
// =============================================================================
Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
GL::Constraint::MinDist::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new GL::Constraint::MinDist(*this, type));
}
// =============================================================================
int
GL::Constraint::MinDist::numConstraints() const
{
  return constraints_.numRows();
}
// =============================================================================
void
GL::Constraint::MinDist::setX(const NOX::Abstract::Vector & y)
{ 
  const NOX::Epetra::MultiVector & yE =
          dynamic_cast<NOX::Epetra::MultiVector&>(*(y.createMultiVector(1)));
  psi_ = glSystem_->getGlKomplex()->real2complex( *(yE.getEpetraMultiVector()(0)) );
  
  isValidConstraints_ = false;

  return;
}
// =============================================================================
void
GL::Constraint::MinDist::setParam( int paramID, double val )
{
  paramsVector_[paramID] = val;
  
  // This would usually be set to FALSE, but in this case, the constraints
  // do not depend of the constraints. Hence, they remain valid.
  // TODO are you sure it must be set to false??
  isValidConstraints_ = false;
}
// =============================================================================
void
GL::Constraint::MinDist::setParams( const vector<int> & paramIDs,
                                const NOX::Abstract::MultiVector::DenseMatrix & vals)
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    paramsVector_[paramIDs[i]] = vals(i,0);

  // This would usually be set to FALSE, but in this case, the constraints
  // do not depend of the constraints. Hence, they remain valid.
  // TODO are you sure it must be set to false??
  isValidConstraints_ = false;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
GL::Constraint::MinDist::computeConstraints()
{
  if (!isValidConstraints_)
  {
    // TODO This goes horribly wrong as Tpetra really calculates
    // Re( psiRef_^T psi_ ) when doing dot().
    // Reported to mailing list on 03/03/2010.
    constraints_(0,0) = std::imag( psiRef_->dot(*psi_) );
    isValidConstraints_ = true;
  }

  return NOX::Abstract::Group::Ok;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
GL::Constraint::MinDist::computeDX()
{
  // TODO compute something here?!
  return NOX::Abstract::Group::Ok;
}
// =============================================================================
// The first column of dgdp should be filled with the constraint
// residuals g if isValidG is false.
// If isValidG is true, then the dgdp contains g on input.
NOX::Abstract::Group::ReturnType
GL::Constraint::MinDist::computeDP( const vector<int> & paramIDs,
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
GL::Constraint::MinDist::isConstraints() const
{
  return isValidConstraints_;
}
// =============================================================================
bool
GL::Constraint::MinDist::isDX() const
{
  return true;
}
// =============================================================================
const NOX::Abstract::MultiVector::DenseMatrix &
GL::Constraint::MinDist::getConstraints() const
{
  return constraints_;
}
// =============================================================================
const NOX::Abstract::MultiVector*
GL::Constraint::MinDist::getDX() const
{
  // TODO check if this is correct
  return NULL;
}
// =============================================================================
bool
GL::Constraint::MinDist::isDXZero() const
{
  // TODO check if this is correct
  return true;
}
// =============================================================================
void
GL::Constraint::MinDist::postProcessContinuationStep( LOCA::Abstract::Iterator::StepStatus stepStatus )
{
  // If the step has been successful, put the solution in psiRef.
  // This means updating the phase condition each *continuation step anew.
  // This is not at each Newton step, which would be rather what we want.
  if (stepStatus == LOCA::Abstract::Iterator::Successful)
    psiRef_ = psi_;

  return;
}
// =============================================================================
