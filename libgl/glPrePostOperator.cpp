#include "glPrePostOperator.h"

#include <Teuchos_ParameterList.hpp>
#include <NOX_Solver_Generic.H>


GlPrePostOperator::GlPrePostOperator(const NOX::Utils& u) :
  numRunPreIterate(0),
  numRunPostIterate(0),
  numRunPreSolve(0),
  numRunPostSolve(0)
{ 
  utils = u;
}

GlPrePostOperator::~GlPrePostOperator()
{
}

void GlPrePostOperator::
runPreIterate(const NOX::Solver::Generic& solver)
{
  ++numRunPreIterate;
  utils.out(NOX::Utils::Details) << 
    "1Dfem's runPreIterate() routine called!" << endl;
}

void GlPrePostOperator::
runPostIterate(const NOX::Solver::Generic& solver)
{
  ++numRunPostIterate;
  utils.out(NOX::Utils::Details) 
    << "1Dfem's runPostIterate() routine called!" << endl;
}

void GlPrePostOperator::
runPreSolve(const NOX::Solver::Generic& solver)
{
  ++numRunPreSolve;
  utils.out(NOX::Utils::Details) 
    << "1Dfem's runPreSolve() routine called!" << endl;
}

void GlPrePostOperator::
runPostSolve(const NOX::Solver::Generic& solver)
{
  ++numRunPostSolve;
  utils.out(NOX::Utils::Details) 
    << "1Dfem's runPostSolve() routine called!" << endl;
}