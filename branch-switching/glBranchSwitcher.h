#ifndef GLBRANCHSWITCHER_H
#define GLBRANCHSWITCHER_H

#include <NOX_Epetra_Interface_Jacobian.H>


#include "glSystem.h"

#include "ginzburgLandau.h"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>

#include <Epetra_CrsMatrix.h>

#include <Teuchos_ParameterList.hpp>

#include <NOX_Epetra_Interface_Required.H> // NOX base class
#include <NOX_Epetra_Interface_Jacobian.H> // NOX base class

#include <LOCA_Epetra_Interface_Required.H> // LOCA base class

#include <NOX_Abstract_PrePostOperator.H>

#include <LOCA_Parameter_Vector.H>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <NOX_Abstract_Group.H>


class GlBranchSwitcher:
      public NOX::Epetra::Interface::Jacobian
  {
  public:

    GlBranchSwitcher ( const Teuchos::RCP<GlSystem> glSystem,
                       const double ds,
                       const double alpha0,
                       const double alpham1,
                       const Teuchos::RCP<const Epetra_Vector> u0,
                       const Teuchos::RCP<const Epetra_Vector> um1 );
                       

    //! Destructor
    ~GlBranchSwitcher();

    //! Evaluate the Ginzburg--Landau functions at a given state defined
    //! by the input vector x.
    virtual bool
    computeF ( const Epetra_Vector &x,
               Epetra_Vector       &F,
               const NOX::Epetra::Interface::Required::FillType fillFlag = NOX::Epetra::Interface::Required::Residual
             );

    //! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
    //! at a given state defined by the input vector x.
    virtual bool
    computeJacobian ( const Epetra_Vector &x,
                      Epetra_Operator     &Jac );

    //! Dummy preconditioner function. So far does nothing but throwing
    //! an exception when called.
    virtual bool
    computePreconditioner ( const Epetra_Vector     &x,
                            Epetra_Operator         &Prec,
                            Teuchos::ParameterList* precParams=0
                          )  const;

    //! Returns the current state. Not necessarily a solution to the problem!
    //! @return Reference-counted pointer to the current state.
    Teuchos::RCP<Epetra_Vector> getSolution() const;

    //! Returns the current Jacobian.
    //! @return Reference-counted pointer to the Jacobian.
    Teuchos::RCP<Epetra_CrsMatrix> getJacobian() const;

  private:

    void
    computePhi2();

    void
    computedFdalpha();

    Teuchos::RCP<Epetra_Vector>
    getNullvector( const Teuchos::RCP<const Epetra_Operator> A );

    const Teuchos::RCP<GlSystem> glSystem_;
    const Teuchos::RCP<const Tpetra::MultiVector<double,int> > phi2_;
    const Teuchos::RCP<const Epetra_Vector > u0_;
    const Teuchos::RCP<const Epetra_Vector > um1_;
    const Teuchos::RCP<const Epetra_Vector>                    dFdalpha_;
    const double ds_;
    const double alpha0_;
    const double alpham1_;

  };

#endif // GLBRANCHSWITCHER_H
