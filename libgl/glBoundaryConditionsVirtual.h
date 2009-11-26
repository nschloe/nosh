#ifndef GLBOUNDARYCONDITIONSVIRTUAL_H
#define GLBOUNDARYCONDITIONSVIRTUAL_H

#include "GridSquare.h"
#include "MagneticVectorPotential.h"

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
    getGlEntry ( const int                                eqIndex,
                 const Tpetra::Vector<double_complex,int> &psi,
                 const GridSquare::GridSquare             &grid,
                 const MagneticVectorPotential            &A
               ) const = 0; // pure virtual

    //! Returns entries and positions of the Jacobian matrix belonging to the
    //! boundary conditions.
    virtual void
    getGlJacobianRow ( const int                                               eqIndex,
                       const Teuchos::RCP<Tpetra::Vector<double_complex,int> > &psi,
                       const GridSquare::GridSquare                            &grid,
                       const MagneticVectorPotential                           &A,
                       const bool                                              fillValues,
                       std::vector<int>                                        &columnIndicesPsi,
                       std::vector<double_complex>                             &valuesPsi,
                       std::vector<int>                                        &columnIndicesPsiConj,
                       std::vector<double_complex>                             &valuesPsiConj
                     ) const = 0; // pure virtual

  protected:

    //! These are very much the grid node types, as appearing the GridVirtual.
    //! The reason for this is that the discretizations applied yield localized
    //! equations which can be associated with exactly one grid point.
    enum equationType
    {
      BOTTOMLEFTCONVEX,
      BOTTOMLEFTCONCAVE,
      BOTTOMRIGHTCONVEX,
      BOTTOMRIGHTCONCAVE,
      TOPLEFTCONVEX,
      TOPLEFTCONCAVE,
      TOPRIGHTCONVEX,
      TOPRIGHTCONCAVE,
      BOTTOM,
      RIGHT,
      TOP,
      LEFT
    };

    //! With the \cboundaryEquationIndex-th boundary equation
    static GlBoundaryConditionsVirtual::equationType
    getEquationType ( const int                      boundaryEquationIndex,
                      const GridSquare::GridSquare & grid
                    );

  private:

  };
#endif // GLBOUNDARYCONDITIONSVIRTUAL_H
