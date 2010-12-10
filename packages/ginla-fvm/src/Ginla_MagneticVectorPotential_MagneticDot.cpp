#include "Ginla_MagneticVectorPotential_MagneticDot.h"

#include "Ginla_EpetraFVM_StkMesh.h"

// ============================================================================
Ginla::MagneticVectorPotential::MagneticDot::
MagneticDot( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
             double mu ) :
  Virtual( mesh ),
  mu_( mu ),
  magnetRadius_( 4.5315 ),
  zz1_ ( 0.1 ),
  zz2_ ( 1.1 ),
  edgeMidpointProjectionCache_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( mesh->getOwnedCells().size() ) ),
  edgeMidpointProjectionCacheUpToDate_( false )
{
}
// =============================================================================
Ginla::MagneticVectorPotential::MagneticDot::
~MagneticDot()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::MagneticDot::
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
Ginla::MagneticVectorPotential::MagneticDot::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "mu", mu_ );

  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::MagneticDot::
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
Ginla::MagneticVectorPotential::MagneticDot::
getRawA_( const Point & x ) const
{
    // Span a cartesian grid over the sample, and integrate over it.

    // The number of grids within the diameter of the magnetic dot, counting
    // the boundary nodes in.
    int n_mag = 101;
    // grid width
    double dx = 2.0 * magnetRadius_ / (n_mag - 1);

    double ax = 0.0;
    double ay = 0.0;

    for ( int ix=0; ix < n_mag; ix++ )
    {
        double xi = ix*dx - magnetRadius_ ;
        for ( int iy=0; iy < n_mag; iy++ )
        {
            double yi = iy*dx - magnetRadius_;
            // circular shape magnetic dot
            if ( xi*xi + yi*yi <= magnetRadius_*magnetRadius_ )
            {
                // x distance between grid point x to magnetic point xi
                double xx = x[0] - xi;
                // y distance between grid point y to magnetic point yi
                double yy = x[1] - yi;
                // r distance between grid point X to magnetic point (xi,yi)
                double r = xx * xx + yy * yy;

                if ( r > 1.0e-15 )
                {
                    // 3D distance to point on upper edge (xi,yi,zz1)
                    double r_3D1 = sqrt( r + zz1_ * zz1_);
                    // 3D distance to point on lower edge (xi,yi,zz2)
                    double r_3D2 = sqrt( r + zz2_ * zz2_);
                    double alpha = ( zz2_ / r_3D2 - zz1_ / r_3D1 ) * dx * dx / r;
                    ax += yy * alpha;
                    ay -= xx * alpha;
                }
            }
        }
    }

    return Teuchos::rcp( new Point( Teuchos::tuple<double>( ax, ay, 0.0 ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::MagneticDot::
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
Ginla::MagneticVectorPotential::MagneticDot::
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
      unsigned int cellDimension = mesh_->getCellDimension( k );
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
