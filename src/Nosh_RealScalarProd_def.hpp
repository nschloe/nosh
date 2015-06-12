#ifndef NOSH_REALSCALARPROD_DEF_HPP
#define NOSH_REALSCALARPROD_DEF_HPP

#include "Nosh_RealScalarProd_decl.hpp"

#include <Thyra_ScalarProdBase.hpp>
#include <Thyra_MultiVectorStdOps.hpp>

namespace Nosh {

template<class Scalar>
bool
RealScalarProd<Scalar>::
isEuclideanImpl() const
{
  return true;
}

template<class Scalar>
void
RealScalarProd<Scalar>::
scalarProdsImpl(
    const Thyra::MultiVectorBase<Scalar>& X,
    const Thyra::MultiVectorBase<Scalar>& Y,
    const Teuchos::ArrayView<Scalar> &scalarProds_out
    ) const
{
  std::cout << "> RealScalarProd<Scalar>::scalarProdsImpl()" << std::endl;
  Thyra::dots(X, Y, scalarProds_out);
  for (int k = 0; k < scalarProds_out.size(); k++) {
    scalarProds_out[k] = getRealPart(scalarProds_out[k]);
  }
  std::cout << "  RealScalarProd<Scalar>::scalarProdsImpl() >" << std::endl;
}

} // end namespace Nosh

#endif  // NOSH_REALSCALARPROD_DEF_HPP
