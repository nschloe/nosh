// @HEADER
//
//    Builds the kinetic energy operator and its variants.
//    Copyright (C) 2010--2012  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
// =============================================================================
// includes
#include "Nosh_MatrixBuilder_Keo.hpp"
#include "Nosh_StkMesh.hpp"
#include "Nosh_ScalarField_Virtual.hpp"
#include "Nosh_VectorField_Virtual.hpp"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>

#include <ml_epetra_preconditioner.h>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
#endif

namespace Nosh {
namespace MatrixBuilder {
// =============================================================================
Keo::
Keo(const Teuchos::RCP<const Nosh::StkMesh> &mesh,
    const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness,
    const Teuchos::RCP<const Nosh::VectorField::Virtual> &mvp
    ) :
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  keoFillTime_( Teuchos::TimeMonitor::getNewTimer(
                  "Nosh: Keo::fillKeo_" ) ),
#endif
  mesh_( mesh ),
  thickness_( thickness ),
  mvp_( mvp ),
  globalIndexCache_(),
  globalIndexCacheUpToDate_( false ),
  keoGraph_(this->buildKeoGraph_()), // build the graph immediately
  keoCache_(Copy, keoGraph_),
  keoBuildParameters_(),
  keoDpCache_(Copy, keoGraph_),
  alphaCache_(),
  alphaCacheUpToDate_( false ),
  paramName_()
{
}
// =============================================================================
Keo::
~Keo()
{
}
// =============================================================================
const Epetra_Comm &
Keo::
getComm() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  return mesh_->getComm();
}
// =============================================================================
const Epetra_FECrsGraph &
Keo::
getGraph() const
{
  return keoGraph_;
}
// =============================================================================
void
Keo::
apply(const std::map<std::string, double> & params,
      const Epetra_Vector &X,
      Epetra_Vector &Y
      ) const
{
  // Check if the relevant parameters from the previous build have changed.
  bool needsRebuild;
  if (keoBuildParameters_.empty())
  {
    needsRebuild = true;
  }
  else
  {
    needsRebuild = false;
    for (std::map<std::string, double>::const_iterator it = keoBuildParameters_.begin();
         it != keoBuildParameters_.end();
         ++it)
    {
      // Check if it->first is in params at all and if their values are equal.
      std::map<std::string, double>::const_iterator it2 = params.find(it->first);
      TEUCHOS_ASSERT(it2 != params.end());
      if (it2->second != it->second)
      {
        needsRebuild = true;
        break;
      }
    }
  }

  // Rebuild if necessary.
  if (needsRebuild)
  {
    this->fillKeo_(keoCache_, params, &Keo::fillerRegular_);
    // Reset build parameters.
    for (std::map<std::string, double>::iterator it = keoBuildParameters_.begin();
         it != keoBuildParameters_.end();
         ++it)
    {
      std::map<std::string, double>::const_iterator it2 = params.find(it->first);
      TEUCHOS_ASSERT(it2 != params.end());
      it->second = it2->second;
    }
  }
  // This direct application in the cache saves
  // one matrix copy compared to an explicit
  // fill(), Apply().
  TEUCHOS_ASSERT_EQUALITY(0, keoCache_.Apply(X, Y));
  return;
}
// =============================================================================
void
Keo::
applyDKDp(const std::map<std::string, double> & params,
          const std::string & paramName,
          const Epetra_Vector &X,
          Epetra_Vector &Y
          ) const
{
  // Unconditionally rebuild the cache.
  // It would be little effort to also wrap this into a conditional
  // with mvpParams and parmamIndex, but it never seems to occur
  // that dK/dp with the same parameter needs to be applied
  // twice in a row.
  paramName_ = paramName;
  this->fillKeo_(keoDpCache_, params, &Keo::fillerDp_);

  TEUCHOS_ASSERT_EQUALITY(0, keoDpCache_.Apply(X, Y));
  return;
}
// =============================================================================
void
Keo::
fill(Epetra_FECrsMatrix &matrix,
     const std::map<std::string, double> & params
     ) const
{
  // Cache the construction of the KEO.
  // This is useful because in the continuation context,
  // getKeo() is called a number of times with the same arguments
  // (in computeF, getJacobian(), and getPreconditioner().
  bool needsRebuild;
  if (keoBuildParameters_.empty())
  {
    needsRebuild = true;
  }
  else
  {
    needsRebuild = false;
    for (std::map<std::string, double>::const_iterator it = keoBuildParameters_.begin();
         it != keoBuildParameters_.end();
         ++it)
    {
      // Check if it->first is in params at all and if their values are equal.
      std::map<std::string, double>::const_iterator it2 = params.find(it->first);
      TEUCHOS_ASSERT(it2 != params.end());
      if (it2->second != it->second)
      {
        needsRebuild = true;
        break;
      }
    }
  }

  if (needsRebuild)
  {
    this->fillKeo_(keoCache_, params, &Keo::fillerRegular_);
    // Reset build parameters.
    for (std::map<std::string, double>::iterator it = keoBuildParameters_.begin();
         it != keoBuildParameters_.end();
         ++it)
    {
      std::map<std::string, double>::const_iterator it2 = params.find(it->first);
      TEUCHOS_ASSERT(it2 != params.end());
      it->second = it2->second;
    }
  }

  matrix = keoCache_;
  return;
}
// =============================================================================
const std::map<std::string,double>
Keo::
getParameters() const
{
  return mvp_->getParameters();
}
// =============================================================================
const Epetra_FECrsGraph
Keo::
buildKeoGraph_() const
{
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the ownedNodes map, and the overlapping one from the
  // overlapNodes.
  // Let's illustrate the implications with the example of the matrix
  //   [ 2 1   ]
  //   [ 1 2 1 ]
  //   [   1 2 ].
  // Suppose subdomain 1 consists of node 1, subdomain 2 of node 3,
  // and node 2 forms the boundary between them.
  // For two processes, if process 1 owns nodes 1 and 2, the matrix
  // will be split as
  //   [ 2 1   ]   [       ]
  //   [ 1 2 1 ] + [       ]
  //   [       ]   [   1 2 ].
  // The vectors always need to have a unique map (otherwise, norms
  // cannot be computed by Epetra), so let's assume they have the
  // map ( [1,2], [3] ).
  // The communucation for a matrix-vector multiplication Ax=y
  // needs to be:
  //
  //   1. Communicate x(3) to process 1.
  //   2. Communicate x(2) to process 2.
  //   3. Compute.
  //
  // If the matrix is split up like
  //   [ 2 1   ]   [       ]
  //   [ 1 1   ] + [   1 1 ]
  //   [       ]   [   1 2 ]
  // (like the overlap map suggests), then any Ax=y comes down to:
  //
  //   1. Communicate x(2) to process 2.
  //   2. Compute.
  //   3. Communicate (part of) y(2) to process 1.
  //
  // In the general case, assuming that the number of nodes adjacent
  // to a boundary (on one side) are approximately the number of
  // nodes on that boundary, there is not much difference in
  // communication between the patterns.
  // What does differ, though, is the workload on the processes
  // during the computation phase: Process 1 that owns the whole
  // boundary, has to compute more than process 2.
  // Notice, however, that the total number of computations is
  // lower in scenario 1 (7 vs. 8 FLOPs); the same is true for
  // storage.
  // Hence, it comes down to the question whether or not the
  // mesh generator provided a fair share of the boundary nodes.
  // If yes, then scenario 1 will yield approximately even
  // computation times; if not, then scenario 2 will guarantee
  // equal computation times at the cost of higher total
  // storage and computation needs.
  //
  // Remark:
  // This matrix will later be fed into ML. ML has certain restrictions as to
  // what maps can be used. One of those is that RowMatrixRowMap() and
  // OperatorRangeMap must be the same, and, if the matrix is square,
  // OperatorRangeMap and OperatorDomainMap must coincide too.
  //
#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  const Epetra_Map &noMap = *mesh_->getComplexNonOverlapMap();
  Epetra_FECrsGraph keoGraph(Copy, noMap, 0);

  const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > edges =
    mesh_->getEdgeNodes();
  if ( !globalIndexCacheUpToDate_ )
    this->buildGlobalIndexCache_( edges );

  // Loop over all edges and put entries wherever two nodes are connected.
  for (unsigned int k=0; k<edges.size(); k++)
    TEUCHOS_ASSERT_EQUALITY( 0,
      keoGraph.InsertGlobalIndices(4, globalIndexCache_[k].Values(),
                                   4, globalIndexCache_[k].Values()));

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  TEUCHOS_ASSERT_EQUALITY(0, keoGraph.GlobalAssemble(noMap,noMap));

  return keoGraph;
}
// =============================================================================
void
Keo::
fillerRegular_(const int k,
               const std::map<std::string, double> & params,
               double * v) const
{
  // Fill v with
  // Re(-exp(i Aint))
  // Im(-exp(i Aint))
  // 1.0
  const double aInt = mvp_->getEdgeProjection(k, params);
  // If compiled with GNU (and maybe other compilers), we could use
  // sincos() here to comute sin and cos simultaneously.
  // PGI, for one, doesn't support sincos, though.
  v[0] = -cos(aInt);
  v[1] = -sin(aInt);
  v[2] = 1.0;
  //sincos( aInt, &sinAInt, &cosAInt );

  return;
}
// =============================================================================
void
Keo::
fillerDp_(const int k,
          const std::map<std::string, double> & params,
          double * v) const
{
  double aInt = mvp_->getEdgeProjection(k, params);
  // paramName_ is set in the KEO building routine.
  double dAdPInt = mvp_->getDEdgeProjectionDp(k, params, paramName_);
  //sincos( aInt, &sinAInt, &cosAInt );
  v[0] =  dAdPInt * sin(aInt);
  v[1] = -dAdPInt * cos(aInt);
  v[2] = 0.0;
  return;
}
// =============================================================================
void
Keo::
fillKeo_(Epetra_FECrsMatrix &keoMatrix,
         const std::map<std::string, double> & params,
         void (Keo::*filler)(const int, const std::map<std::string, double> &, double*) const
         ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm( *keoFillTime_ );
#endif
  // Zero-out the matrix.
  TEUCHOS_ASSERT_EQUALITY(0, keoMatrix.PutScalar(0.0));

#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
#ifndef NDEBUG
  TEUCHOS_ASSERT( !thickness_.is_null() );
  TEUCHOS_ASSERT( !mvp_.is_null() );
#endif

  const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > edges =
    mesh_->getEdgeNodes();
  if ( !globalIndexCacheUpToDate_ )
    this->buildGlobalIndexCache_( edges );
  if ( !alphaCacheUpToDate_ )
    this->buildAlphaCache_(edges, mesh_->getEdgeCoefficients());

  double v[3];
  Epetra_SerialDenseMatrix A(4, 4);
  // Loop over all edges.
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
    // ---------------------------------------------------------------
    // Compute the integral
    //
    //    I = \int_{x0}^{xj} (xj-x0).A(x) / |xj-x0| dx
    //
    // numerically by the midpoint rule, i.e.,
    //
    //    I ~ |xj-x0| * (xj-x0) . A( 0.5*(xj+x0) ) / |xj-x0|.
    //
    // -------------------------------------------------------------------
    // Project vector field onto the edge.
    // Instead of first computing the projection over the normalized edge
    // and then multiply it with the edge length, don't normalize the
    // edge vector.
    (this->*filler)(k, params, v);
    // We'd like to insert the 2x2 matrix
    //
    //     [   alpha                   , - alpha * exp( -IM * aInt ) ]
    //     [ - alpha * exp( IM * aInt ),   alpha                       ]
    //
    // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
    // that shares and edge.
    // Do that now, just blockwise for real and imaginary part.
    v[0] *= alphaCache_[k];
    v[1] *= alphaCache_[k];
    v[2] *= alphaCache_[k];
    A(0, 0) = v[2];   A(0 ,1) = 0.0;     A(0, 2) = v[0];    A(0, 3) = v[1];
    A(1, 0) = 0.0;    A(1 ,1) = v[2];    A(1, 2) = -v[1];   A(1, 3) = v[0];
    A(2, 0) = v[0];   A(2 ,1) = -v[1];   A(2, 2) = v[2];    A(2, 3) = 0.0;
    A(3, 0) = v[1];   A(3 ,1) = v[0];    A(3, 2) = 0.0;     A(3, 3) = v[2];
    TEUCHOS_ASSERT_EQUALITY(0, keoMatrix.SumIntoGlobalValues(
                                                    globalIndexCache_[k], A));
    // -------------------------------------------------------------------
  }

  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix.GlobalAssemble() );
  return;
}
// =============================================================================
void
Keo::
buildGlobalIndexCache_( const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > &edges ) const
{
  globalIndexCache_ =
    Teuchos::ArrayRCP<Epetra_IntSerialDenseVector>(edges.size());

  Teuchos::Tuple<int,2> gid;
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
    gid[0] = edges[k][0]->identifier() - 1;
    gid[1] = edges[k][1]->identifier() - 1;

    globalIndexCache_[k] = Epetra_IntSerialDenseVector( 4 );
    globalIndexCache_[k][0] = 2*gid[0];
    globalIndexCache_[k][1] = 2*gid[0]+1;
    globalIndexCache_[k][2] = 2*gid[1];
    globalIndexCache_[k][3] = 2*gid[1]+1;
  }

  globalIndexCacheUpToDate_ = true;

  return;
}
// =============================================================================
void
Keo::
buildAlphaCache_(const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > & edges,
                 const Teuchos::ArrayRCP<const double> &edgeCoefficients
                 ) const
{
  alphaCache_ = Teuchos::ArrayRCP<double>( edges.size() );

  Teuchos::Tuple<int,2> gid;
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
    // Get the ID of the edge endpoints in the nonoverlapping
    // nodes map.
    gid[0] = edges[k][0]->identifier() - 1;
    gid[1] = edges[k][1]->identifier() - 1;
    const int tlid0 = mesh_->getNodesMap()->LID( gid[0] );
    const int tlid1 = mesh_->getNodesMap()->LID( gid[1] );
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG( tlid0 < 0,
                         "The global index " << gid[0]
                         << " does not seem to be present on this node." );
    TEUCHOS_TEST_FOR_EXCEPT_MSG( tlid1 < 0,
                         "The global index " << gid[1]
                         << " does not seem to be present on this node." );
#endif
    // Update cache.
    std::map<std::string,double> dummy;
    const double thickness = 0.5 * (thickness_->getV(tlid0, dummy) + thickness_->getV(tlid1, dummy));
    alphaCache_[k] = edgeCoefficients[k] * thickness;
  }

  alphaCacheUpToDate_ = true;

  return;
}
// =============================================================================
} // namespace MatrixBuilder
} // namespace Nosh
