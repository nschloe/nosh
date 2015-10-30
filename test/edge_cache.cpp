// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2011--2012  Nico Schl\"omer
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
#include <string>
#include <vector>

#include <Teuchos_ParameterList.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <nosh.hpp>

#include <Teuchos_UnitTestHarness.hpp>

namespace
{

// =============================================================================
void
testCache(const std::string & input_filename_base,
          Teuchos::FancyOStream & out,
          bool & success)
{
  // Read the data from the file.
  auto comm =  Teuchos::DefaultComm<int>::getComm();
  const std::string input_filename = (comm->getSize() == 1) ?
    "data/" + input_filename_base + ".e" :
    "data/" + input_filename_base + "-split.par";

  // Read the data from the file.
  auto mesh = nosh::read(input_filename);

  // Cast the data into something more accessible.
  std::shared_ptr<nosh::Mesh> & mesh = data.get("mesh", std::shared_ptr<nosh::Mesh>());

  std::vectorRCP<double> edge_coefficients;
  std::vectorRCP<DoubleVector> edge_coefficients_fallback;
  unsigned int numEdges;
  if (mesh->supportsEdges()) {
    edge_coefficients = mesh->getEdge_coefficients();
    edge_coefficients_fallback = mesh->getEdge_coefficients_fallback();
    numEdges = mesh->getOverlapEdges().size();
  } else {
    // If we can't get the coefficients, there's nothing to compare.
    return;
  }

  // -------------------------------------------------------------------------
  // Build the equivalent of edge_coefficients using edge_coefficients_fallback.
  std::vectorRCP<double> edge_coefficients2(numEdges);
  std::vector<stk_classic::mesh::Entity*> cells = mesh->getOwnedCells();
  for (unsigned int k = 0; k < cells.size(); k++) {
    stk_classic::mesh::PairIterRelation local_nodes =
      cells[k]->relations(mesh->getMetaData()->node_rank());
    unsigned int num_local_nodes = local_nodes.size();

    stk_classic::mesh::PairIterRelation local_edges =
      cells[k]->relations(mesh->getMetaData()->edge_rank());

    // Fetch the nodal positions into 'local_nodes'.
    const std::vectorRCP<const DoubleVector> local_node_coords =
      mesh->getNode_coordinates(local_nodes);

    // Gather the edge coordinates.
    int num_local_edges = num_local_nodes*(num_local_nodes-1) / 2;
    std::vectorRCP<DoubleVector> local_edge_coords(num_local_edges);
    unsigned int i = 0;
    for (unsigned int e0 = 0; e0 < num_local_nodes; e0++) {
      const int gid0 = (*local_nodes[e0].entity()).identifier() - 1;
      for (unsigned int e1 = e0+1; e1 < num_local_nodes; e1++) {
        const int gid1 = (*local_nodes[e1].entity()).identifier() - 1;

        // Find the edge that has gid0, gid1 as endpoints.
        bool found_edge = false;
        unsigned int j;
        for (j = 0; j < local_edges.size(); j++) {
          // Get the endpoints
          stk_classic::mesh::PairIterRelation end_points =
            (*local_edges[j].entity()).relations(mesh->getMetaData()->node_rank());
          TEUCHOS_ASSERT_EQUALITY(end_points.size(), 2);
          int ep0 = (*end_points[0].entity()).identifier() - 1;
          int ep1 = (*end_points[1].entity()).identifier() - 1;
          if ((ep0 == gid0 && ep1 == gid1) || (ep0 == gid1 && ep1 == gid0)) {
            found_edge = true;
            break;
          }
        }
        TEUCHOS_ASSERT(found_edge);
        // TODO Need LOCAL identifier?
        const int global_edge_id = (*local_edges[j].entity()).identifier() - 1;
        const int local_edge_id = mesh->getEdgesOverlapMap()->LID(global_edge_id);
        TEST_FOR_EXCEPT_MSG(
            local_edge_id < 0,
            "The global index " << global_edge_id
            << " does not seem to be present on this node."
            );

        TEUCHOS_ASSERT_INEQUALITY(edge_coefficients2.size(), >, local_edge_id);
        // Sum the coefficients up.
        edge_coefficients2[local_edge_id] += edge_coefficients_fallback[k][i];
        i++;
      }
    }
  }
  // -------------------------------------------------------------------------
  // Compare edge_coefficients and edge_coefficients2.
  for (unsigned int k = 0; k < numEdges; k++)
    TEST_COMPARE(
      fabs(edge_coefficients[k]-edge_coefficients2[k]),
      <=,
      1.0e-12
    );

//     TEST_COMPARE_FLOATING_ARRAYS(edge_coefficients,
//                                   edge_coefficients2,
//                                   1.0e-12
//                                );
  // -------------------------------------------------------------------------
  return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST(nosh, Edge_cacheRectangleSmall)
{
  testCache("rectanglesmall", out, success);
}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, Edge_cachePacman)
{
  testCache("pacman", out, success);
}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, Edge_cachePacman)
{
  testCache("shell", out, success);
}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, Edge_cachePacman)
{
  testCache("sphere", out, success);
}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, Edge_cache_cubeSmall)
{
  testCache("cubesmall", out, success);
}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, Edge_cacheBrickWHoleHashes)
{
  testCache("brick-w-hole", out, success);
}
// ============================================================================
} // namespace
