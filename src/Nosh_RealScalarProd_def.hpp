#ifndef NOSH_REALSCALARPROD_DEF_HPP
#define NOSH_REALSCALARPROD_DEF_HPP

#include "Nosh_RealScalarProd_decl.hpp"

#include <Thyra_ScalarProdBase.hpp>
#include <Thyra_MultiVectorStdOps.hpp>

namespace nosh {

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
scalar_prodsImpl(
    const Thyra::MultiVectorBase<Scalar>& X,
    const Thyra::MultiVectorBase<Scalar>& Y,
    const Teuchos::ArrayView<Scalar> &scalar_prods_out
    ) const
{
  std::cout << "> RealScalarProd<Scalar>::scalar_prodsImpl()" << std::endl;
  Thyra::dots(X, Y, scalar_prods_out);
  for (int k = 0; k < scalar_prods_out.size(); k++) {
    scalar_prods_out[k] = getRealPart(scalar_prods_out[k]);
  }
  std::cout << "  RealScalarProd<Scalar>::scalar_prodsImpl() >" << std::endl;
}

} // end namespace nosh

#endif  // NOSH_REALSCALARPROD_DEF_HPP
