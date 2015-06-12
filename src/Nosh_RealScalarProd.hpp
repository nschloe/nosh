#include "Nosh_RealScalarProd_decl.hpp"
#include "Nosh_RealScalarProd_def.hpp"

//#include <Thyra_ScalarProdBase.hpp>
//
//namespace Nosh {
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
//        scalarProdsImpl(
//            const Thyra::MultiVectorBase<Scalar> &X,
//            const Thyra::MultiVectorBase<Scalar> &Y,
//            const Teuchos::ArrayView<Scalar> &scalarProds_out
//            ) const
//        {
//          std::cout << "> NOSH::scalarProdsImpl()" << std::endl;
//          Thyra::dots(X, Y, scalarProds_out);
//          //for (int k = 0; k < scalarProds_out.size(); k++) {
//          //  scalarProds_out[k] = getRealPart(scalarProds_out[k]);
//          //}
//          std::cout << "  NOSH::scalarProdsImpl() >" << std::endl;
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
