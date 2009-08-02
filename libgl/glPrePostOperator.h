#include <NOX_Common.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Utils.H>

#include "glSystem.h"

class GlPrePostOperator : public NOX::Abstract::PrePostOperator {

public:

  //! Ctor.
 GlPrePostOperator( Teuchos::RCP<GlSystem>        glsystem,
                    const Teuchos::ParameterList& problemParams);

  //! Destructor.
  ~GlPrePostOperator();

  void runPreIterate(const NOX::Solver::Generic& solver);

protected:

  NOX::Utils utils;

  int numRunPreIterate;

  Teuchos::ParameterList problemParameters_;
  Teuchos::RCP<GlSystem> glsystem_;

};