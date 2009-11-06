#ifndef GLBOUNDARYCONDITIONSCENTRAL_H
#define GLBOUNDARYCONDITIONSCENTRAL_H

#include "glBoundaryConditionsVirtual.h"

#include <complex>
typedef std::complex<double> double_complex;

class GlBoundaryConditionsCentral: public GlBoundaryConditionsVirtual
  {
  public:

    //! Default constructor.
    GlBoundaryConditionsCentral();

    //! Destructor
    virtual ~GlBoundaryConditionsCentral();

    //! Return the value of the Ginzburg-Landau equations for the equation
    //! at eqType.
    virtual double_complex
    getGlEntry ( const int                                eqIndex,
                 const Tpetra::Vector<double_complex,int> &psi,
                 const Grid::Grid                         &grid,
                 const MagneticVectorPotential            &A
               ) const;

    //! Returns entries and positions of the Jacobian matrix belonging to the
    //! boundary conditions.
    virtual void
    getGlJacobianRow ( const int                                               eqIndex,
                       const Teuchos::RCP<Tpetra::Vector<double_complex,int> > &psi,
                       const Grid::Grid                                        &grid,
                       const MagneticVectorPotential                           &A,
                       const bool                                              fillValues,
                       std::vector<int>                                        &columnIndicesPsi,
                       std::vector<double_complex>                             &valuesPsi,
                       std::vector<int>                                        &columnIndicesPsiConj,
                       std::vector<double_complex>                             &valuesPsiConj
                     ) const;

  protected:
  private:

  };
#endif // GLBOUNDARYCONDITIONSCENTRAL_H
