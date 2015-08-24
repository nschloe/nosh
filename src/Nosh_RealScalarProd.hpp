#include "Nosh_RealScalarProd_decl.hpp"
#include "Nosh_RealScalarProd_def.hpp"

//#include <Thyra_ScalarProdBase.hpp>
//
//namespace nosh {
//
//  template<class Scalar>
//    class RealScalarProduct: public Thyra::ScalarProdBase<Scalar>
//  {
//    public:
//      virtual
//        bool
//        isEuclideanImpl() const
//        {
//          return false;
//        };
//
//      virtual
//        void
//        scalar_prodsImpl(
//            const Thyra::MultiVectorBase<Scalar> &X,
//            const Thyra::MultiVectorBase<Scalar> &Y,
//            const Teuchos::ArrayView<Scalar> &scalar_prods_out
//            ) const
//        {
//          std::cout << "> NOSH::scalar_prodsImpl()" << std::endl;
//          Thyra::dots(X, Y, scalar_prods_out);
//          //for (int k = 0; k < scalar_prods_out.size(); k++) {
//          //  scalar_prods_out[k] = getRealPart(scalar_prods_out[k]);
//          //}
//          std::cout << "  NOSH::scalar_prodsImpl() >" << std::endl;
//        }
//
//    private:
//      Scalar
//        getRealPart(Scalar r) const
//        {
//          return r;
//        }
//
//      Scalar getRealPart(std::complex<Scalar> c) const
//      {
//        return c.real();
//      }
//  };
//
//}
