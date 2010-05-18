/*
 * glNoxHelpers.h
 *
 *  Created on: Jan 18, 2010
 *      Author: nico
 */

#ifndef GLNOXHELPERS_H_
#define GLNOXHELPERS_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Comm.h>
#include <NOX_Epetra.H>

#include "Ginla_LocaSystem_Bordered.h"
#include "Recti_Grid_Uniform.h"

namespace glNoxHelpers
{

void
createGlSystem ( const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                 const Teuchos::RCP<const Epetra_Comm>         & eComm,
                 const std::string                             & fileName,
                 Teuchos::ParameterList                        & problemParameters,
                 Teuchos::RCP<Ginla::LocaSystem::Virtual>      & glSystem,
                 Teuchos::RCP<Ginla::State>                    & initialState,
                 Teuchos::RCP<Recti::Grid::Uniform>            & grid
               );

void
createGlSystem ( const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                 const Teuchos::RCP<const Epetra_Comm>         & eComm,
                 const unsigned int Nx,
                 const double scaling,
                 const double H0,
                 const Teuchos::ParameterList               & domainParameters,
                 Teuchos::ParameterList                     & problemParameters,
                 Teuchos::RCP<Ginla::LocaSystem::Virtual>   & glSystem,
                 Teuchos::RCP<Ginla::State>                 & initialState,
                 Teuchos::RCP<Recti::Grid::Uniform>         & grid
               );

Teuchos::RCP<NOX::Epetra::Group>
createSolverGroup ( const Teuchos::RCP<Ginla::LocaSystem::Virtual> & glSystem,
                    const Teuchos::RCP<Teuchos::ParameterList>     & nlParamsPtr,
                    const Teuchos::RCP<Ginla::State>               & initialState
                  );

Teuchos::RCP<NOX::StatusTest::Generic>
createConvergenceTest ( Teuchos::ParameterList & noxStatusList,
                        Teuchos::ParameterList & nlParamsPtr );

Teuchos::RCP<NOX::Solver::Generic>
createSolver ( const Teuchos::RCP<NOX::Epetra::Group>       grpPtr,
               const Teuchos::RCP<NOX::StatusTest::Generic> statusTest,
               const Teuchos::RCP<Teuchos::ParameterList>   nlParamsPtr );

double
computeJacobianConditionNumber ( const Teuchos::RCP<const NOX::Solver::Generic> solver,
                                 const Teuchos::RCP<      NOX::Epetra::Group>   grpPtr );

void
computeJacobianEigenvalues ( const Teuchos::RCP<const NOX::Solver::Generic> solver,
                             const Teuchos::RCP<      NOX::Epetra::Group>   grpPtr,
                             const int MyPID );

void
printSolutionToFile ( const std::string                                    & outputDir,
                      const std::string                                    & fileBaseName,
                      const std::string                                    & outputFormat,
                      const Teuchos::RCP<const NOX::Solver::Generic>       & solver,
                      const Teuchos::RCP<const Ginla::LocaSystem::Virtual> & glSystem
                    );

int
checkConvergence ( const Teuchos::RCP<const NOX::Solver::Generic> solver );

}

#endif /* GLNOXHELPERS_H_ */
