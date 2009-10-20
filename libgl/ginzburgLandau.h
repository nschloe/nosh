/********************************************//**
 * The Ginzburg--Landau equations.
 ***********************************************/
#ifndef GINZBURGLANDAU_H
#define GINZBURGLANDAU_H

#include <complex>
#include <string> // for the file output functions

#include "glBoundaryConditionsVirtual.h"
#include "staggeredGrid.h"

#include <Teuchos_ParameterList.hpp>

#include <Teuchos_XMLObject.hpp>

#include <Tpetra_MultiVector.hpp>

#include <Epetra_Map.h>

// abbreviate the complex type name
typedef std::complex<double> double_complex;

class GinzburgLandau
{
  public:

     /*! Default constructor.
         @param nx Number of boxes in each spatial dimension.
         @param edgelength Edge length of the square domain.
         @param h0 (Initial) external magnetic field strength. */
     GinzburgLandau( int nx,
                     double edgelength,
                     double h0,
                     Teuchos::RCP<GlBoundaryConditionsVirtual> bc );

     /*! Destructor. */
     ~GinzburgLandau();

     /*! Returns a pointer to the \f$A\f$ grid in use.*/
     StaggeredGrid::StaggeredGrid* getStaggeredGrid();
     
     Tpetra::MultiVector<double_complex,int>
     computeGlVector( Tpetra::MultiVector<double_complex,int> psi );

     /*! Returns the coefficients of the jacobian system associated with the
         Ginzburg--Landau equations. */
     void
     getJacobianRow( const int                                     eqnum,
                     const Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi,
                     std::vector<int>                              &columnIndicesPsi,
                     std::vector<double_complex>                   &valuesPsi,
                     std::vector<int>                              &columnIndicesPsiConj,
		     std::vector<double_complex>                   &valuesPsiConj );

     /*! Get sparsity pattern of the jacobian system. */
     void
     getJacobianRowSparsity( int              eqnum,
                             std::vector<int> &columnIndicesPsi,
                             std::vector<int> &columnIndicesPsiConj );

     /*! Calcuate the grid approximation of the Gibbs free energy
       \f[
       \mathcal{G} = \int\nolimits_{\Omega} |\psi|^4 \,\mathrm{d}\omega
       \f]
       of a given state \f$\psi\f$. */
     double
     freeEnergy ( const Tpetra::MultiVector<double_complex,int> &psi );
     
     /*! Calculate the vorticity of the current solution. */
     int
     getVorticity ( const Tpetra::MultiVector<double_complex,int> &psi );
      
  private:

     /*! Evaluates the Ginzburg--Landau equations.
         @param eqnum Index of the equation to evaluate.
         @param psi   Current order parameter \f$\psi\f$. */
     double_complex
     computeGl( const int                                     eqnum,
                const Tpetra::MultiVector<double_complex,int> &psi   );
    
      //! Equation type enumerator.
      /*! Semantically separates the different types of conditions which must
          be applied at different parts of the rectangular grid. */
      enum equationType
      {
        BOUNDARY,
        INTERIOR,
        PHASE_CONDITION
      };

      StaggeredGrid::StaggeredGrid sGrid;

      void getEquationType ( const int           eqnum,
                             equationType        &eqType,
                             int                 &eqIndex );
      
      /*! Have \c boundaryConditions_ declared as pointer as its class
          \c GlBoundaryConditionsVirtual is only available via a forward
	  declaration. */
      Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions_;

      /*! Calculated the coefficients of the jacobian system associated with the
          Ginzburg--Landau equations. */
      void computeJacobianRow ( const bool                                    fillValues,
                                const int                                     eqnum,
                                const Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi,
                                std::vector<int>                              &columnIndicesPsi,
                                std::vector<double_complex>                   &valuesPsi,
                                std::vector<int>                              &columnIndicesPsiConj,
                                std::vector<double_complex>                   &valuesPsiConj );

};
#endif // GINZBURGLANDAU_H