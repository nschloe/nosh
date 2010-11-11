#ifndef GINLA_FVM_KINETICENERGYOPERATOR_H
#define GINLA_FVM_KINETICENERGYOPERATOR_H
// =============================================================================
#include <Tpetra_Operator.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "VIO_TpetraMesh_Mesh.h"
#include "Ginla_MagneticVectorPotential_Virtual.h"
#include "Ginla_Typedefs.h"
// =============================================================================
//typedef MultiVector<double_complex> ComplexMultivector;
typedef Tpetra::CrsGraph<ORD> TCrsGraph;
typedef Tpetra::CrsMatrix<double_complex,ORD> ComplexMatrix;
// =============================================================================
namespace Ginla {
namespace FVM {

class KineticEnergyOperator : public Tpetra::Operator
{
public:
    KineticEnergyOperator( const Teuchos::RCP<VIO::TpetraMesh::Mesh>                   & mesh,
                           const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp
                         );

    virtual
    ~KineticEnergyOperator();

    //! Returns the Map associated with the domain of this operator,
    //! which must be compatible with X.getMap().
    virtual
    const Teuchos::RCP<const Map> &
    getDomainMap () const;

    //! Returns the Map associated with the range of this operator,
    //! which must be compatible with Y.getMap().
    virtual
    const Teuchos::RCP <Map> &
    getRangeMap () const ;

    //! Computes the operator-multivector application.
    virtual
    void apply ( const ComplexMultiVector & X,
                       ComplexMultiVector & Y,
                 Teuchos::ETransp mode=Teuchos::NO_TRANS,
                 double_complex alpha = 1.0,
                 double_complex beta  = 0.0
               ) const;


    //! Indicates whether this operator supports applying the adjoint operator.
    virtual
    bool hasTransposeApply () const;

protected:
private:

    void
    assembleKeo_() const;

    void
    createKeoGraph_() const;

private:

//    bool useTranspose_;
//    const Teuchos::RCP<const Epetra_Comm> comm_;

    const Teuchos::RCP<VIO::TpetraMesh::Mesh> mesh_;
    const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp_;

    mutable Teuchos::RCP<Tpetra::CrsGraph> keoGraph_;
    mutable Teuchos::RCP<Tpetra::CrsMatrix> keo_;

//    mutable double keoMu_;
//    mutable Teuchos::Tuple<double,3> keoScaling_;
//
//    Teuchos::RCP<Epetra_LinearProblem> keoProblem_;
//    mutable Teuchos::RCP<Amesos_BaseSolver> keoSolver_;

};

} // namespace FVM
} // namespace Ginla
// =============================================================================
#endif // GINLA_FVM_KINETICENERGYOPERATOR_H
