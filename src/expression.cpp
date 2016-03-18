#include "expression.hpp"

namespace nosh {

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  integrate_over_control_volumes(
      const nosh::expression & expr,
      const nosh::mesh & mesh
      )
  {
    auto cv = mesh.control_volumes();
    auto vals = std::make_shared<Tpetra::Vector<double,int,int>>(cv->getMap());

    // can only integrate degree 0 for now
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        expr.degree > 0,
        "Expression degree " << expr.degree << " too high. "
        << "Can only integrate expressions of degree 0 for now."
        );

    auto cv_data = cv->getData();
    auto vals_data = vals->getDataNonConst();
    auto nodes = mesh.get_owned_nodes();

#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(nodes.size(), cv_data.size());
    TEUCHOS_ASSERT_EQUALITY(nodes.size(), vals_data.size());
#endif

    switch (expr.degree) {
      case 0:
        for (size_t k = 0; k < nodes.size(); k++)
        {
          vals_data[k] = expr(mesh.get_coords(nodes[k])) * cv_data[k];
        }
        break;

      default:
        throw "Illegal degree";
    }

    return vals;
  }

} // namespace nosh
