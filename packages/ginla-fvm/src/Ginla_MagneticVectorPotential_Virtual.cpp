#include "Ginla_MagneticVectorPotential_Virtual.h"

#include "Ginla_EpetraFVM_StkMesh.h"

// ============================================================================
Ginla::MagneticVectorPotential::Virtual::
Virtual( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh ):
  mesh_( mesh )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::Virtual::
~Virtual()
{
}
// ============================================================================
