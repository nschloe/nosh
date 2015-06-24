#include "Expression.hpp"

namespace Nosh {

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  integrateOverControlVolumes(
      const Nosh::Expression & expr,
      const Nosh::Mesh & mesh
      )
  {
    auto cv = mesh.getControlVolumes();
    auto vals = std::make_shared<Tpetra::Vector<double,int,int>>(cv->getMap());

    // can only integrate degree 0 for now
    TEUCHOS_ASSERT_INEQUALITY(expr.degree, <, 1);

    auto cvData = cv->getData();
    auto valsData = vals->getDataNonConst();
    const VectorFieldType & coordsField = mesh.getNodeField("coordinates");
    auto nodes = mesh.getOwnedNodes();

#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(nodes.size(), cvData.size());
    TEUCHOS_ASSERT_EQUALITY(nodes.size(), valsData.size());
#endif

    switch (expr.degree) {
      case 0:
        for (size_t k = 0; k < nodes.size(); k++)
        {
          auto x = mesh.getNodeValue(coordsField, nodes[k]);
          valsData[k] = expr.eval(x) * cvData[k];
        }
        break;

      default:
        throw "Illegal degree";
    }

    return vals;
  }

} // namespace Nosh
