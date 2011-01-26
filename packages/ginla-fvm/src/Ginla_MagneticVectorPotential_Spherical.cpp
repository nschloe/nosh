#include "Ginla_MagneticVectorPotential_Spherical.h"

#include "Ginla_EpetraFVM_StkMesh.h"

// ============================================================================
Ginla::MagneticVectorPotential::Spherical::
Spherical( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
           double mu,
           double phi,
           double theta ):
  Virtual( mesh ),
  mu_( mu ),
  phi_( phi ),
  theta_( theta ),
  edgeMidpointProjectionCache_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( mesh->getOwnedCells().size() ) ),
  edgeMidpointProjectionCacheUpToDate_( false )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::Spherical::
~Spherical()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::Spherical::
setParameters( const LOCA::ParameterVector & p )
{
    bool valuesChanged = false;

    if (p.isParameter( "mu" ))
        if ( mu_ != p.getValue ( "mu" ) )
        {
            mu_ = p.getValue ( "mu" );
            valuesChanged = true;
        }

    if (p.isParameter( "phi" ))
        if ( phi_ != p.getValue ( "phi" ) )
        {
            phi_ = p.getValue ( "phi" );
            valuesChanged = true;
        }

    if (p.isParameter( "theta" ))
        if ( theta_ != p.getValue ( "theta" ) )
        {
            theta_ = p.getValue ( "theta" );
            valuesChanged = true;
        }

    return valuesChanged;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::MagneticVectorPotential::Spherical::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "mu", mu_ );
  p->addParameter( "phi", phi_ );
  p->addParameter( "theta", theta_ );

  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Spherical::
getA(const Point & x) const
{
  Teuchos::RCP<Point> a = this->getRawA_( x );
  (*a)[0] *= mu_;
  (*a)[1] *= mu_;
  (*a)[2] *= mu_;
  return a;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Spherical::
getRawA_(const Point & x) const
{
  return
  Teuchos::rcp( new Point( Teuchos::tuple<double>(
    -0.5 * sin(theta_)             * x[1] + 0.5 * cos(theta_) * sin(phi_) * x[2],
     0.5 * sin(theta_)             * x[0] - 0.5 * cos(theta_) * cos(phi_) * x[2],
     0.5 * cos(theta_) * cos(phi_) * x[1] - 0.5 * cos(theta_) * sin(phi_) * x[0]
  ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Spherical::
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
Ginla::MagneticVectorPotential::Spherical::
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

              // Get the midpoint of the edge.
              Point midpoint; // get A(midpoint)
              for (int i=0; i<midpoint.size(); i++ )
                  midpoint[i] = 0.5 * ( node0[i] + node1[i] );

              Teuchos::RCP<Point> a = this->getRawA_( midpoint );

              // Instead of first computing the projection over the normalized edge
              // and then multiply it with the edge length, don't normalize the
              // edge vector.
              edgeMidpointProjectionCache_[k][edgeIndex] = 0.0;
              for (int i=0; i<midpoint.size(); i++ )
                  edgeMidpointProjectionCache_[k][edgeIndex] += ( node1[i] - node0[i] ) * (*a)[i];
              edgeIndex++;
          }
      }
  }

  edgeMidpointProjectionCacheUpToDate_ = true;
  return;
}
// ============================================================================