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

#include "GlMinDistConstraint.h"

#include "glSystem.h"

// =============================================================================
// Default constructor
GlMinDistConstraint::GlMinDistConstraint ( const Teuchos::RCP <GlSystem> & glSystem,
                                           const NOX::Abstract::Vector & initialGuess ):
  glSystem_( glSystem ),
  constraints_(1,1),
  isValidConstraints_(false),
  paramsVector_(glSystem->GetContinuableParams()),
  psi_(Teuchos::null),
  psiRef_(Teuchos::null)
{ 
  // Initialize constraint
  constraints_.putScalar(0.0);

  // At the first step, the reference solution is
  // the initial guess. For other steps, the solution
  // at the previous continuation step will be taken 
  // (see postProcessContinuationStep())
  psiRef_ = glSystem_->getGlKomplex()->real2complex( initialGuess );
  
  // Soluion equals initialGuess at the instantiation.
  // Essentially we needed to shape x as xRef
  psi_ = glSystem_->getGlKomplex()->real2complex( initialGuess );
}
// =============================================================================
// Destructor
GlMinDistConstraint::~GlMinDistConstraint()
{
}
// =============================================================================
// copy constructor
GlMinDistConstraint::GlMinDistConstraint(const GlMinDistConstraint& source, 
                                         NOX::CopyType              type) :
  glSystem_(source.glSystem_),
  constraints_(source.constraints_),
  isValidConstraints_(false),
  psi_(source.psi_->clone(type)),
  psiRef_(source.psiRef_->clone(type))
{
  if (source.isValidConstraints_ && type == NOX::DeepCopy)
    isValidConstraints_ = true;
}
// =============================================================================
GlMinDistConstraint::~GlMinDistConstraint()
{
}
// =============================================================================
void
GlMinDistConstraint::copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const GlMinDistConstraint& source = dynamic_cast<const GlMinDistConstraint&>(src);

  if (this != &source) {
    glSystem_ = source.glSystem_;
    constraints_ = source.constraints_;
    isValidConstraints_ = source.isValidConstraints_;
    *psi_ = *source.psi_;
    *psiRef_ = *source.psiRef_;
  }

  return;
}
// =============================================================================
Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
GlMinDistConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new GlMinDistConstraint(*this, type));
}
// =============================================================================
int
GlMinDistConstraint::numConstraints() const
{
  return constraints_.numRows();
}
// =============================================================================
void
GlMinDistConstraint::setX(const NOX::Abstract::Vector & y)
{
  psi_ = glSystem_->getGlKomplex()->real2complex( y );
  isValidConstraints_ = false;

  return;
}
// =============================================================================
void
GlMinDistConstraint::setParam( int paramID, double val )
{
  paramsVector_[paramID] = val;
  
  // This would usually be set to FALSE, but in this case, the constraints
  // do not depend of the constraints. Hence, they remain valid.
  // isValidConstraints_ = false;
}
// =============================================================================
void
GlMinDistConstraint::setParams( const vector<int> & paramIDs, 
                                const NOX::Abstract::MultiVector::DenseMatrix & vals)
{
  for (int i=0; i<paramIDs.size(); i++)
    paramsVector_[paramIDs[i]] = vals(i,0);

  // This would usually be set to FALSE, but in this case, the constraints
  // do not depend of the constraints. Hence, they remain valid.
  // isValidConstraints_ = false;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
GlMinDistConstraint::computeConstraints()
{
  if (!isValidConstraints_)
  {
    // TODO Check if this incorporates Hermitian conjugation!
    constraints_(0,0) = std::imag( psiRef_->dot(psi_) );
    isValidConstraints_ = true;
  }

  return NOX::Abstract::Group::Ok;
}
// =============================================================================
NOX::Abstract::Group::ReturnType
GlMinDistConstraint::computeDX()
{
  // TODO compute something here?!
  return NOX::Abstract::Group::Ok;
}
// =============================================================================
// The first column of dgdp should be filled with the constraint
// residuals g if isValidG is false.
// If isValidG is true, then the dgdp contains g on input.
NOX::Abstract::Group::ReturnType
GlMinDistConstraint::computeDP( const vector<int> & paramIDs, 
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
GlMinDistConstraint::isConstraints() const
{
  return isValidConstraints_;
}
// =============================================================================
bool
GlMinDistConstraint::isDX() const
{
  return true;
}
// =============================================================================
const NOX::Abstract::MultiVector::DenseMatrix &
GlMinDistConstraint::getConstraints() const
{
  return constraints_;
}
// =============================================================================
const NOX::Abstract::MultiVector*
GlMinDistConstraint::getDX() const
{
  // TODO check if this is correct
  return NULL;
}
// =============================================================================
bool
GlMinDistConstraint::isDXZero() const
{
  // TODO check if this is correct
  return true;
}
// =============================================================================
void
GlMinDistConstraint::postProcessContinuationStep( LOCA::Abstract::Iterator::StepStatus stepStatus )
{
  // If the step has been successful, put the solution in psiRef.
  // This means updating the phase condition each *continuation step anew.
  // This is not at each Newton step, which would be rather what we want.
  if (stepStatus == LOCA::Abstract::Iterator::Successful)
    psiRef_ = psi_;

  return;
}
// =============================================================================