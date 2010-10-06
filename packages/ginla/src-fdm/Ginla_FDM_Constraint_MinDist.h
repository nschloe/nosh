/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2009--2010 Nico Schl\"omer

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

//! This class implements the \c LOCA::MultiContinuation::ConstraintInterfaceMVDX
//! constraint interface for the phase constraint
//! \f[
//! g(\psi) = \Im (\tilde{\psi}^{\mathrm{H}} \psi)
//! \f]
//! where \f$\psi\f$ is the solution vector, and \f$\tilde{\psi}\f$ a
//! reference solution.
//! This constraint does not depend upon the constrained parameter \f$\mu\f$.

#ifndef GLMINDISTCONSTRAINT
#define GLMINDISTCONSTRAINT

#include <LOCA_MultiContinuation_ConstraintInterface.H>

#include "Ginla_FDM_LocaSystem_Default.h"

#include <LOCA_Parameter_Vector.H>

namespace Ginla {
  
  namespace FDM {

  namespace Constraint {

class MinDist:
            // Can't use LOCA::MultiContinuation::ConstraintInterfaceMVDX
            // as the multiplication with the constraint is calculated
            // incorrectly, i.e., there's no way to define a vector
            // such that multiplication with it equals our phase condition.
            public LOCA::MultiContinuation::ConstraintInterface
{
  public:

  // Constructor
  MinDist( const Teuchos::RCP<const Ginla::Komplex::LinearProblem> & komplex,
           const Teuchos::RCP<ComplexVector>                       & psi,
           const LOCA::ParameterVector                             & paramsVector
         );

  // Copy constructor
  MinDist( const Ginla::FDM::Constraint::MinDist & source,
                 NOX::CopyType                     type = NOX::DeepCopy
         );

  // Destructor
  virtual
  ~MinDist();

  // Copy
  virtual void 
  copy(const LOCA::MultiContinuation::ConstraintInterface& source);

  // Cloning function
  virtual Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
  clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Before the next continuation step, save the solution vector.
  virtual void
  preProcessContinuationStep( LOCA::Abstract::Iterator::StepStatus stepStatus );
  
  // Return number of constraints
  virtual int
  numConstraints() const;

  // Set the solution vector to y.
  virtual void
  setX(const NOX::Abstract::Vector& y);

  // Sets parameter indexed by paramID
  virtual void
  setParam(int paramID, double val);

  // Sets parameters indexed by paramIDs
  virtual void
  setParams( const vector<int>& paramIDs, 
             const NOX::Abstract::MultiVector::DenseMatrix& vals);

  // Compute continuation constraint equations
  virtual NOX::Abstract::Group::ReturnType
  computeConstraints();

  // Compute derivative of constraints w.r.t. solution vector x
  virtual NOX::Abstract::Group::ReturnType
  computeDX();

  //! Compute derivative of constraints w.r.t. supplied parameters.
  virtual NOX::Abstract::Group::ReturnType
  computeDP(const vector<int>& paramIDGlMinDistConstraints, 
            NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
            bool isValidG);

  // Return true if constraint residuals are valid
  virtual bool
  isConstraints() const;
  
  // Return true if derivative of constraint w.r.t. x is valid
  virtual bool
  isDX() const;

  // Return constraint residuals
  virtual const NOX::Abstract::MultiVector::DenseMatrix&
  getConstraints() const;

  //! Return \c true if solution component of constraint 
  //! derivatives is zero.
  virtual bool
  isDXZero() const;
  
  //! Compute result_p = alpha * dg/dx * input_x.
  virtual NOX::Abstract::Group::ReturnType
  multiplyDX ( double                                    alpha,
               const NOX::Abstract::MultiVector        & input_x,
               NOX::Abstract::MultiVector::DenseMatrix & result_p
             ) const;
             
  //! Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x.
  virtual NOX::Abstract::Group::ReturnType        
  addDX ( Teuchos::ETransp transb,
          double alpha,
          const NOX::Abstract::MultiVector::DenseMatrix &b,
          double beta,
          NOX::Abstract::MultiVector &result_x
        ) const;
  
private:

  // Prohibit generation and use of operator=()
  MinDist& operator=(const MinDist& source);

private:
  
  //! ProblemLOCAPrototype problem
  Teuchos::RCP<const Ginla::Komplex::LinearProblem> komplex_;

  //! Constraint values
  NOX::Abstract::MultiVector::DenseMatrix constraints_;

  //! Flag indicating whether constraints are valid (up-to-date)
  bool isValidConstraints_;

  //! LOCA parameter vector
  LOCA::ParameterVector paramsVector_;

  //! Solution vector
  Teuchos::RCP<ComplexVector> psi_;

  //! Reference Solution vector
  const Teuchos::RCP<ComplexVector> psiRef_;
};

    } // namespace Constraint
  }
} // namespace GL
#endif // GLMINDISTCONSTRAINT