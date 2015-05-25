#include "Ginla_FVM_KineticEnergyOperator.h"

// =============================================================================
Ginla::FVM::KineticEnergyOperator::
KineticEnergyOperator( const Teuchos::RCP<VIO::TpetraMesh::Mesh>                   & mesh,
                       const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp
                     ):
    mesh_ ( mesh ),
    mvp_( mvp ),
    keo_( Teuchos::null ),
    keoGraph_( Teuchos::null )
{
}
// =============================================================================
Ginla::FVM::KineticEnergyOperator::
~KineticEnergyOperator()
{
}
// =============================================================================
const Teuchos::RCP<const Map> &
Ginla::FVM::KineticEnergyOperator::
getDomainMap () const
{
    TEUCHOS_ASSERT( !keo_.is_null() );
    return keo_->getDomainMap();
}
// =============================================================================
const Teuchos::RCP <Map> &
Ginla::FVM::KineticEnergyOperator::
getRangeMap () const
{
    return keo_->getRangeMap();
}
// =============================================================================
void
Ginla::FVM::KineticEnergyOperator::
apply ( const ComplexMultiVector & X,
              ComplexMultiVector & Y,
        Teuchos::ETransp mode,
        double_complex alpha,
        double_complex beta
      ) const
{
    if ( keo_.is_null() )
        this->assembleKeo_();

    return keo_->apply( X, Y, mode, alpha, beta );
}
// =============================================================================
bool
Ginla::FVM::KineticEnergyOperator::
hasTransposeApply () const
{
    // operator is Hermitian anyway
    return true;
}
// =============================================================================
void
Ginla::FVM::KineticEnergyOperator::
assembleKeo_() const
{
      // TODO Don't throw away the old matrix.
      // This whole things could be treated more elegantly with an FEParameterMatrix class.
      // Chris Baker, Aug 22, 2010:
      // In Tpetra, this could happen via some non-member function or helper class.
      // For example, FEParameterMatrix would accept your assembled matrix and internally
      // stored a matrix of non-locals. All of the methods to FEParameterMatrix would send
      // local entries to your assembled matrix and cache non-local entries in another matrix.
      // When you call FEParameterMatrix::doneBuilding(), then it would call the import/export
      // functionality. The graph is only rebuilt if necessary (e.g., because the set of
      // occupied columns grows).
      int maxNumEntriesPerRow = komplex_->getComplexMap()->getGlobalNumElements();

      keo_ = Teuchos::rcp( new ComplexMatrix( komplex_->getComplexMap(),
                                              maxNumEntriesPerRow
                                            )
                         );
      kineticEnergyOperator_->setAllToScalar( 0.0 );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // Loop over the elements, create local load vector and mass matrix,
      // and insert them into the global matrix.
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems = mesh_->getElems();
      Teuchos::ArrayRCP<Point> nodes = mesh_->getNodes();
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeLengths = mesh_->getEdgeLengths();
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > coedgeLengths = mesh_->getCoedgeLengths();

      TEUCHOS_ASSERT( !elems.is_null() );
      TEUCHOS_ASSERT( !nodes.is_null() );
      TEUCHOS_ASSERT( !edgeLengths.is_null() );
      TEUCHOS_ASSERT( !coedgeLengths.is_null() );

      // loop over the local elements
      for ( int k=0; k<elems.size(); k++ )
      {
          Teuchos::ArrayRCP<ORD> & elem = elems[k];

          // loop over the edges
          int n = elem.size();
          for ( int l=0; l<n; l++ )
          {
              // two subsequent indices determine an edge
              Teuchos::Tuple<ORD,2> nodeIndices;
              nodeIndices[0] = elem[l];
              nodeIndices[1] = elem[(l+1)%n];

              // co-edge / edge ratio
              double alpha = coedgeLengths[k][l]
                           / edgeLengths[k][l];

              // -------------------------------------------------------------------
              // Compute the integral
              //
              //    I = \int_{x0}^{xj} (xj-x0).A(x) dx
              //
              // numerically by the midpoint rule, i.e.,
              //
              //    I ~ |xj-x0| * (xj-x0) . A( 0.5*(xj+x0) ).
              //
              Point midpoint; // get A(midpoint)
              Point & node0 = nodes[ nodeIndices[0] ];
              Point & node1 = nodes[ nodeIndices[1] ];
              for (int i=0; i<midpoint.size(); i++ )
                  midpoint[i] = 0.5 * ( node0[i] + node1[i] );

              Teuchos::Tuple<double_complex,2> values;
              // -------------------------------------------------------------------
              // Project vector field onto the edge.
              Teuchos::RCP<Point> a = mvp_->getA( midpoint );
              double aProj = 0.0;
              for (int i=0; i<midpoint.size(); i++ )
                  aProj += ( node1[i] - node0[i] ) * (*a)[i];

              const double aInt = aProj * edgeLengths[k][l];

              values[0] =   alpha;
              values[1] = - alpha * exp( -IM * aInt );
              kineticEnergyOperator_->insertGlobalValues( nodeIndices[0],
                                                          nodeIndices,
                                                          values
                                                        );

              values[0] = - alpha * exp(  IM * aInt ); // integration the other way around
              values[1] =   alpha;
              kineticEnergyOperator_->insertGlobalValues( nodeIndices[1],
                                                          nodeIndices,
                                                          values
                                                        );
          }
      }

      kineticEnergyOperator_->globalAssemble();
      kineticEnergyOperator_->fillComplete();

      return;
}
// =============================================================================
void
Ginla::FVM::KineticEnergyOperator::
createKeoGraph_() const
{
  keoGraph_ = Teuchos::rcp( new TCrsGraph( mesh_->getComplexValuesMap()) );

  TEUCHOS_ASSERT( !mesh_.is_null() );
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems = mesh_->getElems();

  // connect all nodes with all other nodes in each element
  for ( int k=0; k<elems.size(); k++ )
  {
      Teuchos::ArrayRCP<ORD> & elem = elems[k];
      for ( int l=0; l<elem.size(); l++ )
          kineticEnergyOperatorGraph_->insertGlobalIndices( elem[l], elem() );
  }

  keoGraph_->FillComplete();

  return;
}
// =============================================================================
