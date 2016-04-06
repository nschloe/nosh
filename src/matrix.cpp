#include "matrix.hpp"

namespace nosh {
// apply boundary conditions to the matrix
void
matrix::
apply_bcs_()
{
  const auto boundary_vertices = this->mesh->boundary_vertices();
  for (const auto boundary_vertex: boundary_vertices) {
    // check if any of the boundary conditions kicks in
    const auto coord = this->mesh->get_coords(boundary_vertex);
    int count = 0;
    for (const auto & bc: bcs) {
      if (bc->is_inside(coord)) {
        count++;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        count > 1,
        "More than one active boundary conditions. Abort."
        );
    if (count == 1) {
      // eliminate the row in A
      const auto gid = this->mesh->gid(boundary_vertex);
      size_t num = this->getNumEntriesInGlobalRow(gid);
      // It shouldn't actually happen that the specified global row
      // does not belong to this graph.
      // TODO find out if why we need this, fix the underlying issue, make this
      // a TEUCHOS_TEST_*
      if (num != Teuchos::OrdinalTraits<size_t>::invalid()) {
        std::vector<int> cols(num);
        std::vector<double> vals(num);
        this->getGlobalRowCopy(gid, cols, vals, num);
        // set vals to 0
        std::fill(vals.begin(), vals.end(), 0.0);
        // set diagonal entry to 1
        auto it = std::find(cols.begin(), cols.end(), gid);
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            it == cols.end(),
            "Matrix has no main diagonal entry."
            );
        int pos = it - cols.begin();
        vals[pos] = 1.0;
        this->replaceGlobalValues(gid, cols, vals);
      }
    }
  }
}

} // namespace nosh
