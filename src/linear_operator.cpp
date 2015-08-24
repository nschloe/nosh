
#include "linear_operator.hpp"
namespace nosh{

// apply boundary conditions to the matrix
void
linear_operator::
apply_bcs_()
{
  const auto boundary_nodes = this->mesh->boundary_nodes();
  const vector_fieldType & coords_field = this->mesh->get_node_field("coordinates");
  for (const auto boundary_node: boundary_nodes) {
    // check if any of the boundary conditions kicks in
    const auto coord = this->mesh->get_node_value(coords_field, boundary_node);
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
      const auto gid = this->mesh->gid(boundary_node);
      size_t num = this->getNumEntriesInGlobalRow(gid);
      std::vector<int> cols(num);
      std::vector<double> vals(num);
      this->getGlobalRowCopy(gid, cols, vals, num);
      // set vals to 9
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

} // namespace nosh
