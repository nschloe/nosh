#include <NOX_Common.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Utils.H>

#include "glSystem.h"

class GlPrePostOperator : public NOX::Abstract::PrePostOperator {

public:

  //! Constructor.
 GlPrePostOperator( Teuchos::RCP<GlSystem>        glsystem,
                    const Teuchos::ParameterList& problemParams);

  //! Destructor.
  ~GlPrePostOperator();

  //! Function that gets called before each iteration.
  //! This particular implementation prints the current state to the file
  //! data/newton-step-numRunPreIterate.vtk .
  //! @param solver The solver.
  void runPreIterate(const NOX::Solver::Generic& solver);

protected:

  //! How ofter the function has been invoked yet.
  int numRunPreIterate;

  Teuchos::ParameterList problemParameters_;   //!< The problem parameters.
  Teuchos::RCP<GlSystem> glsystem_;            //!< The Ginzburg--Landau system.

};