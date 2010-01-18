/*
 * glNoxHelpers.h
 *
 *  Created on: Jan 18, 2010
 *      Author: nico
 */

#ifndef GLNOXHELPERS_H_
#define GLNOXHELPERS_H_

#include <Teuchos_RCP.hpp>
//#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Comm.h>

//#include <Tpetra_Vector.hpp>
//
//#include <NOX.H>
//#include <NOX_Epetra.H>
//
//#include <NOX.H>
//#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_Epetra.H>
//
//// for the eigenvalue computation:
//#include <AnasaziBasicEigenproblem.hpp>
//#include <AnasaziBasicOutputManager.hpp>
//#include <AnasaziBlockDavidsonSolMgr.hpp>
//#include <AnasaziBlockKrylovSchurSolMgr.hpp>
//#include <AnasaziLOBPCGSolMgr.hpp>
//#include <AnasaziRTRSolMgr.hpp>
//#include <AnasaziConfigDefs.hpp>
//#include <AnasaziEpetraAdapter.hpp>
//#include <AnasaziBasicSort.hpp>
//
//#include <Teuchos_XMLParameterListHelpers.hpp>
//
//
//#include "glNox.h"
//
//#include "ioVirtual.h"
//#include "ioFactory.h"
//
//#include "glBoundaryConditionsVirtual.h"
//#include "glBoundaryConditionsInner.h"
//#include "glBoundaryConditionsOuter.h"
//#include "glBoundaryConditionsCentral.h"
//#include "ginzburgLandau.h"
//#include "glPrePostOperator.h"
//#include "GridUniformSquare.h"
//
//#include "GridSquare.h"
#include "GlSystemWithConstraint.h"

namespace glNoxHelpers {
// =========================================================================
void
createGlSystem( const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
		        const Teuchos::RCP<const Epetra_Comm>         & eComm,
		        const std::string                             & fileName,
		              Teuchos::ParameterList                  & problemParameters,
		              Teuchos::RCP<GlSystemWithConstraint>    & glSystem );

void
createGlSystem( const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                const Teuchos::RCP<const Epetra_Comm>         & eComm,
                const unsigned int Nx,
                const double scaling,
                const double H0,
	                  Teuchos::ParameterList                  & problemParameters,
                      Teuchos::RCP<GlSystemWithConstraint>    & glSystem );

void
setPrePostWriter(       Teuchos::ParameterList & noxParaList,
		          const Teuchos::RCP<const AbstractStateWriter> & asw,
                  const std::string            & outputDir );

Teuchos::RCP<NOX::Epetra::Group>
createSolverGroup( const Teuchos::RCP<GlSystemWithConstraint> glSystem,
		           const Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr );

Teuchos::RCP<NOX::StatusTest::Generic>
createConvergenceTest( Teuchos::ParameterList & noxStatusList,
		               Teuchos::ParameterList & nlParamsPtr );

Teuchos::RCP<NOX::Solver::Generic>
createSolver( const Teuchos::RCP<NOX::Epetra::Group>       grpPtr,
		      const Teuchos::RCP<NOX::StatusTest::Generic> statusTest,
		      const Teuchos::RCP<Teuchos::ParameterList>   nlParamsPtr );

double
computeJacobianConditionNumber( const Teuchos::RCP<const NOX::Solver::Generic> solver,
		                        const Teuchos::RCP<      NOX::Epetra::Group>   grpPtr );

void
computeJacobianEigenvalues( const Teuchos::RCP<const NOX::Solver::Generic> solver,
		                    const Teuchos::RCP<      NOX::Epetra::Group>   grpPtr,
		                    const int MyPID );

void
printSolutionToFile( const Teuchos::RCP<const NOX::Solver::Generic> solver,
		             const Teuchos::RCP<const GlSystemWithConstraint> glSystem,
		             const std::string & fileName );

int
checkConvergence( const Teuchos::RCP<const NOX::Solver::Generic> solver );
}
#endif /* GLNOXHELPERS_H_ */
