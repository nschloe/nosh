#include "Ginla_MagneticVectorPotential_Virtual.h"

// ============================================================================
Ginla::MagneticVectorPotential::Virtual::
Virtual( double mu ) :
  mu_( mu )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::Virtual::
~Virtual()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::Virtual::
setMu( const double mu )
{
  bool valuesChanged = false;

  if ( mu_ != mu )
  {
      mu_ = mu;
      valuesChanged = true;
  }

  return valuesChanged;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Virtual::
getMu() const
{
    return mu_;
}
// ============================================================================
