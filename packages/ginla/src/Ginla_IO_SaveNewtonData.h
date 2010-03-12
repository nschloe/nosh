#ifndef GLPREPOSTOPERATOR_H
#define GLPREPOSTOPERATOR_H

#include <NOX_Common.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>

#include "Ginla_IO_StateWriter.h"
#include "Ginla_Komplex.h"

namespace Ginla {
  namespace IO {

class SaveNewtonData:
        public NOX::Abstract::PrePostOperator
{

public:

  //! Constructor.
  SaveNewtonData ( const Teuchos::RCP<const Ginla::IO::StateWriter> & stateWriter,
                   const Teuchos::RCP<const Recti::Grid::General>   & grid,
                   const Teuchos::RCP<const Ginla::Komplex>         & komplex
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

  const Teuchos::RCP<const Ginla::IO::StateWriter> stateWriter_;
  const Teuchos::RCP<const Recti::Grid::General>   grid_;
  const Teuchos::RCP<const Ginla::Komplex>         komplex_;
};

  } // namespace IO
} // namespace SaveNewtonData

#endif // GLPREPOSTOPERATOR_H
