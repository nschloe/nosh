#include "Ginla_MagneticVectorPotential_Custom.h"

#include "Ginla_EpetraFVM_StkMesh.h"
#include <Epetra_Vector.h>

// ============================================================================
Ginla::MagneticVectorPotential::Custom::
Custom( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
        const Teuchos::RCP<const Epetra_MultiVector>  & mvp,
        double mu
      ):
  Virtual( mesh ),
  mvp_( mvp ),
  mu_( mu ),
  edgeMidpointProjectionCache_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( mesh->getOwnedCells().size() ) ),
  edgeMidpointProjectionCacheUpToDate_( false )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::Custom::
~Custom()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::Custom::
setParameters( const LOCA::ParameterVector & p )
{
    bool valuesChanged = false;

    if (p.isParameter( "mu" ))
        if ( mu_ != p.getValue ( "mu" ) )
        {
            mu_ = p.getValue ( "mu" );
            valuesChanged = true;
        }

    return valuesChanged;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::MagneticVectorPotential::Custom::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "mu", mu_ );

  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Custom::
getA(const Point & x) const
{
  TEST_FOR_EXCEPTION( true,
                      std::runtime_error,
                      "getA(x) for general 'x' cannot be implemented as no general MVP info is provided to this class."
                    );
  Teuchos::RCP<Point> a;;
  return a;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Custom::
getAEdgeMidpointProjection( const unsigned int cellIndex,
                            const unsigned int edgeIndex
                          ) const
{
  if ( !edgeMidpointProjectionCacheUpToDate_ )
      this->initializeEdgeMidpointProjectionCache_();

  return mu_ * edgeMidpointProjectionCache_[cellIndex][edgeIndex];
}
// ============================================================================
void
Ginla::MagneticVectorPotential::Custom::
initializeEdgeMidpointProjectionCache_() const
{
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges and create the cache.
  // To this end, loop over all cells and the edges within the cell.
  for ( int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();

      unsigned int numLocalNodes = rel.size();
      unsigned int cellDimension = mesh_->getCellDimension( numLocalNodes );
      // extract the nodal coordinates
      Teuchos::Array<Point> localNodes = mesh_->getNodeCoordinates( rel );

      edgeMidpointProjectionCache_[k] = Teuchos::ArrayRCP<double>( mesh_->getNumEdgesPerCell( cellDimension ) );

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      Teuchos::Tuple<int,2> gid;
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
          const Point & node0 = localNodes[e0];
          gid[0] = (*rel[e0].entity()).identifier() - 1;
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
              const Point & node1 = localNodes[e1];
              gid[1] = (*rel[e1].entity()).identifier() - 1;

              // Approximate the value at the midpoint of the edge
              // by the average of the values at the adjacent nodes.
              Point a;
              for (int i=0; i<a.size(); i++ )
                  a[i] = 0.5 * ( (*(*mvp_)(i))[gid[0]] + (*(*mvp_)(i))[gid[1]] );

              // Instead of first computing the projection over the normalized edge
              // and then multiply it with the edge length, don't normalize the
              // edge vector.
              edgeMidpointProjectionCache_[k][edgeIndex] = 0.0;
              for (int i=0; i<a.size(); i++ )
                  edgeMidpointProjectionCache_[k][edgeIndex] += ( node1[i] - node0[i] ) * a[i];

              edgeIndex++;
          }
      }
  }

  edgeMidpointProjectionCacheUpToDate_ = true;
  return;
}
// ============================================================================