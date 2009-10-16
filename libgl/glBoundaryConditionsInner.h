#ifndef GLBOUNDARYCONDITIONSINNER_H
#define GLBOUNDARYCONDITIONSINNER_H

#include "glBoundaryConditionsVirtual.h"

#include <Tpetra_MultiVector.hpp>

#include <complex>
typedef std::complex<double> double_complex;

class GlBoundaryConditionsInner: public GlBoundaryConditionsVirtual
  {
  public:

    //! Default constructor.
    GlBoundaryConditionsInner();

    //! Destructor
    virtual ~GlBoundaryConditionsInner();

    //! Return the value of the Ginzburg-Landau equations for the equation
    //! at eqType.
    virtual double_complex
    getGlEntry ( const int                                     eqIndex,
                 const Tpetra::MultiVector<double_complex,int> &psi,
                 const StaggeredGrid::StaggeredGrid            &sGrid
               );

    //! Returns entries and positions of the Jacobian matrix belonging to the
    //! boundary conditions.
    virtual void
    getGlJacobianRow ( const int                                                    eqIndex,
                       const Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi,
                       const StaggeredGrid::StaggeredGrid                           &sGrid,
                       const bool                                                   fillValues,
                       std::vector<int>                                             &columnIndicesPsi,
                       std::vector<double_complex>                                  &valuesPsi,
                       std::vector<int>                                             &columnIndicesPsiConj,
                       std::vector<double_complex>                                  &valuesPsiConj
                     );

  protected:
  private:

  };
#endif // GLBOUNDARYCONDITIONSINNER_H
