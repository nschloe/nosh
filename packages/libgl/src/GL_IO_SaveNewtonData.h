#ifndef GLPREPOSTOPERATOR_H
#define GLPREPOSTOPERATOR_H

#include <NOX_Common.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>

#include "AbstractStateWriter.h"

namespace GL {
  namespace IO {

class SaveNewtonData:
        public NOX::Abstract::PrePostOperator
{

public:

  //! Constructor.
  SaveNewtonData( const Teuchos::RCP<const AbstractStateWriter> & stateWriter,
                  const std::string                             & outputDir 
                );

  //! Destructor.
  ~SaveNewtonData();

  //! Function that gets called before each iteration.
  //! This particular implementation prints the current state to the file
  //! data/newton-step-numRunPreIterate.vtk .
  //! @param solver The solver.
  void runPostIterate(const NOX::Solver::Generic& solver);

protected:
private:
  //! How ofter the function has been invoked yet.
  int numRunPreIterate;

  const Teuchos::RCP<const AbstractStateWriter> stateWriter_;
  std::string                                   outputDir_;
};

  } // namespace IO
} // namespace SaveNewtonData

#endif // GLPREPOSTOPERATOR_H
