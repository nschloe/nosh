/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2009--2010 Nico Schl\"omer, Daniele Avitabile

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

#ifndef GLPERTURBATIONVIRTUAL
#define GLPERTURBATIONVIRTUAL

#include<LOCA_Parameter_Vector.H>
#include<complex>

namespace Ginla {

  namespace Perturbation {

class Virtual
{
  public:

  // Constructor
  Virtual();

  // Destructor
  virtual
  ~Virtual();

  // Compute perturbation
  virtual std::complex<double> 
    computePerturbation( int k ) const = 0;

  // Set parameters
  virtual void 
    setParameters( const LOCA::ParameterVector & p ) = 0;

  protected:

  LOCA::ParameterVector paramList_;

  private:

};

  } // namespace Constraint
} // namespace GL
#endif // GLPERTURBATIONVIRTUAL
