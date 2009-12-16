#ifndef GLPREPOSTOPERATOR_H
#define GLPREPOSTOPERATOR_H

#include <NOX_Common.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>


#include "AbstractStateWriter.h"

class GlPrePostOperator : public NOX::Abstract::PrePostOperator {

public:

  //! Constructor.
  GlPrePostOperator( Teuchos::RCP<AbstractStateWriter> & stateWriter,
                     const Teuchos::ParameterList      & problemParams,
                     const std::string                 & outputDir );

  //! Destructor.
  ~GlPrePostOperator();

  //! Function that gets called before each iteration.
  //! This particular implementation prints the current state to the file
  //! data/newton-step-numRunPreIterate.vtk .
  //! @param solver The solver.
  void runPostIterate(const NOX::Solver::Generic& solver);

protected:

  //! How ofter the function has been invoked yet.
  int numRunPreIterate;

  Teuchos::ParameterList            problemParameters_;   //!< The problem parameters.
  Teuchos::RCP<AbstractStateWriter> stateWriter_;
  std::string                       outputDir_;

};
#endif // GLPREPOSTOPERATOR_H
