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
    const vector_fieldType & coords_field = mesh.get_node_field("coordinates");
    auto nodes = mesh.owned_nodes();

#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(nodes.size(), cv_data.size());
    TEUCHOS_ASSERT_EQUALITY(nodes.size(), vals_data.size());
#endif

    switch (expr.degree) {
      case 0:
        for (size_t k = 0; k < nodes.size(); k++)
        {
          auto x = mesh.get_node_value(coords_field, nodes[k]);
          vals_data[k] = expr.eval(x) * cv_data[k];
        }
        break;

      default:
        throw "Illegal degree";
    }

    return vals;
  }

} // namespace nosh
