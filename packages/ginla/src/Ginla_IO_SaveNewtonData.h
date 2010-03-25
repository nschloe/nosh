#ifndef GLPREPOSTOPERATOR_H
#define GLPREPOSTOPERATOR_H

#include <Teuchos_RCP.hpp>
#include <NOX_Abstract_PrePostOperator.H>

// foward declarations
namespace Ginla {
  class Komplex;
  namespace IO {
      class StateWriter;
  }
}
namespace Recti{
  namespace Grid {
    class General;
  }
}

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
