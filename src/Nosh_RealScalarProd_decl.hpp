#ifndef NOSH_REALSCALARPROD_DECL_HPP
#define NOSH_REALSCALARPROD_DECL_HPP

#include "Thyra_ScalarProdBase_decl.hpp"

namespace nosh {

template<class Scalar>
class RealScalarProd : public Thyra::ScalarProdBase<Scalar> {
protected:
  virtual bool isEuclideanImpl() const;

  virtual void scalar_prodsImpl(
    const Thyra::MultiVectorBase<Scalar>& X,
    const Thyra::MultiVectorBase<Scalar>& Y,
    const Teuchos::ArrayView<Scalar> &scalar_prods
    ) const;

private:

  Scalar
  getRealPart(Scalar r) const
  {
    return r;
  };

  Scalar
  getRealPart(std::complex<Scalar> c) const
  {
    return c.real();
  };

};
} // end namespace nosh

#endif  // NOSH_REALSCALARPROD_DECL_HPP

