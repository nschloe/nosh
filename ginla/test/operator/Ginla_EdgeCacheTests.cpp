// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2011  Nico Schl\"omer
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
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Vector.h>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_StkMesh.hpp"
#include "Ginla_StkMeshReader.hpp"
#include "Ginla_MagneticVectorPotential.hpp"
#include "Ginla_KeoFactory.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// =============================================================================
void
testCache( const std::string & inputFileNameBase,
           Teuchos::FancyOStream & out,
           bool & success )
{
    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
            Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm> eComm =
            Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

    std::string inputFileName;
    if ( eComm->NumProc() == 1 )
        inputFileName = inputFileNameBase + ".e";
    else
        inputFileName = inputFileNameBase + "-balanced.par";
    // =========================================================================
    // Read the data from the file.
    Teuchos::ParameterList data;
    Ginla::StkMeshRead( *eComm, inputFileName, data );

    // Cast the data into something more accessible.
    Teuchos::RCP<Ginla::StkMesh> & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );

    Teuchos::ArrayRCP<double> edgeCoefficients;
    Teuchos::ArrayRCP<DoubleVector> edgeCoefficientsFallback;
    int numEdges;
    try
    {
        edgeCoefficients = mesh->getEdgeCoefficients();
        edgeCoefficientsFallback = mesh->getEdgeCoefficientsFallback();
        numEdges = mesh->getOverlapEdges().size();
    }
    catch( ... )
    {
        // If we can't get the coefficients, there's nothing to compare.
        return;
    }

    // -------------------------------------------------------------------------
    // Build the equivalent of edgeCoefficients using edgeCoefficientsFallback.
    Teuchos::ArrayRCP<double> edgeCoefficients2( numEdges );
    std::vector<stk::mesh::Entity*> cells = mesh->getOwnedCells();
    for (unsigned int k=0; k < cells.size(); k++)
    {
        stk::mesh::PairIterRelation localNodes =
            cells[k]->relations( mesh->getMetaData()->node_rank() );
        unsigned int numLocalNodes = localNodes.size();

        stk::mesh::PairIterRelation localEdges =
            cells[k]->relations( mesh->getMetaData()->edge_rank() );

        // Fetch the nodal positions into 'localNodes'.
        const Teuchos::ArrayRCP<const DoubleVector> localNodeCoords =
            mesh->getNodeCoordinates( localNodes );

        // Gather the edge coordinates.
        int numLocalEdges = numLocalNodes*(numLocalNodes-1) / 2;
        Teuchos::ArrayRCP<DoubleVector> localEdgeCoords(numLocalEdges);
        unsigned int i = 0;
        for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
        {
            const int gid0 = (*localNodes[e0].entity()).identifier() - 1;
            for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
            {
                const int gid1 = (*localNodes[e1].entity()).identifier() - 1;

                // Find the edge that has gid0, gid1 as endpoints.
                bool edgeFound = false;
                int j;
                for ( j=0; j<localEdges.size(); j++ )
                {
                    // Get the endpoints
                    stk::mesh::PairIterRelation endPoints =
                        (*localEdges[j].entity()).relations( mesh->getMetaData()->node_rank() );
                    TEUCHOS_ASSERT_EQUALITY( endPoints.size(), 2 );
                    int ep0 = (*endPoints[0].entity()).identifier() - 1;
                    int ep1 = (*endPoints[1].entity()).identifier() - 1;
                    if ( (ep0==gid0 && ep1==gid1) || (ep0==gid1 && ep1==gid0) )
                    {
                        edgeFound = true;
                        break;
                    }    
                }
                TEUCHOS_ASSERT( edgeFound );
                // TODO Need LOCAL identifier?
                const int gEdgeId = (*localEdges[j].entity()).identifier() - 1;
                const int lEdgeId = mesh->getEdgesOverlapMap()->LID( gEdgeId );
                TEST_FOR_EXCEPT_MSG( lEdgeId < 0,
                                    "The global index " << gEdgeId
                                    << " does not seem to be present on this node." );

//                 if ( gEdge==0 )
//                 {
//                   std::cout << "cell # " << k << std::endl;
//                   std::cout << "local edge # " << i << std::endl;
//                   std::cout << "global edge # " << gEdge << std::endl;
//                   std::cout << "AAA " << edgeCoefficientsFallback[k][i] << std::endl;
//                 }
                TEUCHOS_ASSERT_INEQUALITY( edgeCoefficients2.size(), >, lEdgeId );
                // Sum the coefficients up.
                edgeCoefficients2[lEdgeId] += edgeCoefficientsFallback[k][i];
                i++;
            }
        }
    }
    // -------------------------------------------------------------------------
    // Compare edgeCoefficients and edgeCoefficients2.
    for ( unsigned int k=0; k<numEdges; k++ )
        TEST_COMPARE( fabs(edgeCoefficients[k]-edgeCoefficients2[k]), <=, 1.0e-12 );

//     TEST_COMPARE_FLOATING_ARRAYS( edgeCoefficients,
//                                   edgeCoefficients2,
//                                   1.0e-12
//                                 );
    // -------------------------------------------------------------------------
    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Ginla, EdgeCacheRectangleSmall )
{
    std::string inputFileNameBase = "rectanglesmall";

    testCache( inputFileNameBase,
               out,
               success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, EdgeCachePacman )
{
    std::string inputFileNameBase = "pacman";

    testCache( inputFileNameBase,
               out,
               success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, EdgeCacheCubeSmall )
{
    std::string inputFileNameBase = "cubesmall";

    testCache( inputFileNameBase,
               out,
               success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, EdgeCacheCubeLargeHashes )
{
    std::string inputFileNameBase = "cubelarge";

    testCache( inputFileNameBase,
               out,
               success );
}
// ============================================================================
} // namespace