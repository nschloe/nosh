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

#include "GlSystemLocaConstraint.h"

#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>
#include <NOX_Utils.H>

#include <EpetraExt_RowMatrixOut.h>

#include <EpetraExt_Utils.h>

#include <Epetra_Map.h>

#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Thyra_EpetraThyraWrappers.hpp>

#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// abbreviate the complex type name
typedef std::complex<double> double_complex;

// =============================================================================
// Default constructor
GlSystemLocaConstraint::GlSystemLocaConstraint ( GinzburgLandau::GinzburgLandau &gl,
                                                 const Teuchos::RCP<const Epetra_Comm> eComm,
                                                 const Teuchos::RCP<const ComplexVector> psi,
                                                 const std::string outputDir,
                                                 const std::string outputDataFileName,
                                                 const std::string outputFileFormat,
                                                 const std::string solutionFileNameBase,
                                                 const std::string nullvectorFileNameBase,
                                                 const unsigned int maxStepNumberDecimals
                                               ) :
        glSystem_ ( gl, eComm, psi, outputDir, outputDataFileName, outputFileFormat,
                    solutionFileNameBase, nullvectorFileNameBase, maxStepNumberDecimals )
{
}
// =============================================================================
// Destructor
GlSystemLocaConstraint::~GlSystemLocaConstraint()
{
}
// =============================================================================
bool
GlSystemLocaConstraint::computeF ( const Epetra_Vector & x,
                                   Epetra_Vector       & FVec,
                                   const NOX::Epetra::Interface::Required::FillType fillFlag )
{
    return glSystem_.computeF ( x, FVec, fillFlag );
}
// =============================================================================
bool
GlSystemLocaConstraint::computeJacobian ( const Epetra_Vector & x,
                                          Epetra_Operator     & Jac
                                        )
{
    return glSystem_.computeJacobian ( x, Jac );
}
// =============================================================================
bool
GlSystemLocaConstraint::computePreconditioner ( const Epetra_Vector    & x,
                                                Epetra_Operator        & Prec,
                                                Teuchos::ParameterList * precParams )
{
    return glSystem_.computePreconditioner( x, Prec, precParams );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
GlSystemLocaConstraint::getSolution() const
{
    return glSystem_.getSolution();
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GlSystemLocaConstraint::getJacobian() const
{
    return glSystem_.getJacobian();
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GlSystemLocaConstraint::getPreconditioner() const
{
    return glSystem_.getPreconditioner();
}
// =============================================================================
bool
GlSystemLocaConstraint::computeShiftedMatrix ( double alpha,
                                               double beta,
                                               const Epetra_Vector   & x,
                                                     Epetra_Operator & A )
{   
    return glSystem_.computeShiftedMatrix ( alpha, beta, x , A );
}
// =============================================================================
// function used by LOCA
void
GlSystemLocaConstraint::setParameters ( const LOCA::ParameterVector &p )
{
    glSystem_.setParameters ( p );
}
// =============================================================================
void
GlSystemLocaConstraint::setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper )
{
    glSystem_.setLocaStepper ( stepper );
}
// =============================================================================
void
GlSystemLocaConstraint::releaseLocaStepper()
{
    glSystem_.releaseLocaStepper();
}
// =============================================================================
// function used by LOCA
void
GlSystemLocaConstraint::printSolution ( const  Epetra_Vector &x,
                                        double conParam )
{
    glSystem_.printSolution ( x, conParam );
}
// =============================================================================
// function used by LOCA
void GlSystemLocaConstraint::setOutputDir ( const string &directory )
{
    glSystem_.setOutputDir ( directory );
}
// =============================================================================
void
GlSystemLocaConstraint::writeSolutionToFile ( const Epetra_Vector & x,
                                              const std::string   & filePath
                                            ) const
{
    glSystem_.writeSolutionToFile ( x, filePath );
}
// =============================================================================
void
GlSystemLocaConstraint::writeAbstractStateToFile ( const Epetra_Vector & x,
                                                   const std::string   & filePath
                                                 ) const
{
    glSystem_.writeAbstractStateToFile ( x, filePath );
}
// =============================================================================
// TODO delete?
const Teuchos::RCP<const GlKomplex>
GlSystemLocaConstraint::getGlKomplex() const
{
    return glSystem_.getGlKomplex();
}
// =============================================================================
double
GlSystemLocaConstraint::getH0() const
{
    return glSystem_.getH0();
}
// =============================================================================
void
GlSystemLocaConstraint::setH0 ( const double h0 )
{
    glSystem_.setH0 ( h0 );
    return;
}
// =============================================================================
void
GlSystemLocaConstraint::setScaling ( const double scaling )
{
    glSystem_.setScaling ( scaling );
    return;
}
// =============================================================================
void
GlSystemLocaConstraint::setChi ( const double chi )
{
    glSystem_.setChi ( chi );
    return;
}
// =============================================================================