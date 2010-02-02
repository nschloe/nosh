#ifndef GLBOUNDARYCONDITIONSVIRTUAL_H
#define GLBOUNDARYCONDITIONSVIRTUAL_H

#include "GridUniformVirtual.h"
#include "MagneticVectorPotential.h"

#include <Tpetra_MultiVector.hpp>

#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

#include <complex>
typedef std::complex<double> double_complex;
typedef Tpetra::Vector<double_complex,Thyra::Ordinal> ComplexVector;

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
    getGlEntry ( const int                        eqIndex,
                 const ComplexVector            & psi,
                 const double                     chi,
                 const GridUniformVirtual       & grid,
                 const MagneticVectorPotential  & A
               ) const = 0; // pure virtual

    //! Returns entries and positions of the Jacobian matrix belonging to the
    //! boundary conditions.
    virtual void
    getGlJacobianRow ( const int                           eqIndex,
                       const Teuchos::RCP<ComplexVector> & psi,
                       const double                        chi,
                       const GridUniformVirtual          & grid,
                       const MagneticVectorPotential     & A,
                       const bool                          fillValues,
                       std::vector<int>                  & columnIndicesPsi,
                       std::vector<double_complex>       & valuesPsi,
                       std::vector<int>                  & columnIndicesPsiConj,
                       std::vector<double_complex>       & valuesPsiConj
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
    getEquationType ( const int                  boundaryEquationIndex,
                      const GridUniformVirtual & grid
                    );

  private:

  };
#endif // GLBOUNDARYCONDITIONSVIRTUAL_H
