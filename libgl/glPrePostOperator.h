#include <NOX_Common.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Utils.H>

class GlPrePostOperator : public NOX::Abstract::PrePostOperator {

public:

  //! Ctor.
  GlPrePostOperator(const NOX::Utils& u);

  //! Destructor.
  ~GlPrePostOperator();

  void runPreIterate(const NOX::Solver::Generic& solver);

  void runPostIterate(const NOX::Solver::Generic& solver);

  void runPreSolve(const NOX::Solver::Generic& solver);

  void runPostSolve(const NOX::Solver::Generic& solver);

  int getNumRunPreIterate() const { return numRunPreIterate; };

  int getNumRunPostIterate() const { return numRunPostIterate; };

  int getNumRunPreSolve() const { return numRunPreSolve; };

  int getNumRunPostSolve() const { return numRunPostSolve; };

protected:

  NOX::Utils utils;

  int numRunPreIterate;
  int numRunPostIterate;
  int numRunPreSolve;
  int numRunPostSolve;

};