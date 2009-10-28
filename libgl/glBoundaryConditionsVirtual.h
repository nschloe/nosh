#ifndef GLBOUNDARYCONDITIONSVIRTUAL_H
#define GLBOUNDARYCONDITIONSVIRTUAL_H

#include "staggeredGrid.h"

#include <Tpetra_MultiVector.hpp>

#include <complex>
typedef std::complex<double> double_complex;

class GlBoundaryConditionsVirtual
  {
  public:

    //! Default constructor.
    GlBoundaryConditionsVirtual();

    //! Destructor
    virtual ~GlBoundaryConditionsVirtual();

    //! Return the value of the Ginzburg-Landau equations for the equation
    //! at eqType.
    virtual double_complex
    getGlEntry ( const int                                     eqIndex,
                 const Tpetra::Vector<double_complex,int> &psi,
                 const StaggeredGrid::StaggeredGrid            &sGrid
               ) = 0; // pure virtual

    //! Returns entries and positions of the Jacobian matrix belonging to the
    //! boundary conditions.
    virtual void
    getGlJacobianRow ( const int                                                    eqIndex,
                       const Teuchos::RCP<Tpetra::Vector<double_complex,int> > psi,
                       const StaggeredGrid::StaggeredGrid                           &sGrid,
                       const bool                                                   fillValues,
                       std::vector<int>                                             &columnIndicesPsi,
                       std::vector<double_complex>                                  &valuesPsi,
                       std::vector<int>                                             &columnIndicesPsiConj,
                       std::vector<double_complex>                                  &valuesPsiConj
                     ) = 0; // pure virtual

  protected:

    enum equationType
    {
      BOTTOMLEFT,
      BOTTOMRIGHT,
      TOPRIGHT,
      TOPLEFT,
      BOTTOM,
      RIGHT,
      TOP,
      LEFT
    };

    // TODO: Document.
    static void
    getEquationType ( const int                                 eqIndex,
                      const StaggeredGrid::StaggeredGrid        &sGrid,
                      GlBoundaryConditionsVirtual::equationType &eqType,
                      Teuchos::Array<int>                       &i
                    );

  private:

  };
#endif // GLBOUNDARYCONDITIONSVIRTUAL_H
