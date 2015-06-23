
#include "LinearOperator.hpp"
namespace Nosh{

// apply boundary conditions to the matrix
void
LinearOperator::
applyBcs_()
{
  const auto boundaryNodes = this->mesh->getBoundaryNodes();
  const VectorFieldType & coordsField = this->mesh->getNodeField("coordinates");
  for (const auto boundaryNode: boundaryNodes) {
    // check if any of the boundary conditions kicks in
    const auto coord = this->mesh->getNodeValue(coordsField, boundaryNode);
    int count = 0;
    for (const auto & bc: bcs) {
      if (bc->isInside(coord)) {
        count++;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        count > 1,
        "More than one active boundary conditions. Abort."
        );
    if (count == 1) {
      // eliminate the row in A
      const auto gid = this->mesh->gid(boundaryNode);
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

} // namespace Nosh
