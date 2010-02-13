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
  GlPrePostOperator( const Teuchos::RCP<const AbstractStateWriter> & stateWriter,
                     const std::string                             & outputDir,
                     const std::string                             & outputFormat );

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

  const Teuchos::RCP<const AbstractStateWriter> stateWriter_;
  std::string                                   outputDir_;
  std::string                                   filenameExtension_;

};
#endif // GLPREPOSTOPERATOR_H
