/********************************************//**
 * The Ginzburg--Landau equations.
 ***********************************************/
#include <complex>

#include "psiGrid.h"
#include "aGrid.h"

// abbreviate the complex type name
typedef std::complex<double> double_complex;

class GinzburgLandau
{
  public:

      //! Equation type enumerator.
      /*! Semantically separates the different types of conditions which must
          be applied at different parts of the rectangular grid. */
      enum equationType
      {
        BOTTOMLEFT,
        BOTTOMRIGHT,
        TOPLEFT,
        TOPRIGHT,
        BOTTOM,
        TOP,
        LEFT,
        RIGHT,
        INTERIOR,
        PHASE_CONDITION
      };

     /*! Default constructor. */
     GinzburgLandau( int nx,
                     double edgelength,
                     double h0 );

     /*! Destructor. */
     ~GinzburgLandau();

     /*! Evaluates the Ginzburg--Landau equations. */
     double_complex computeGl( const int                         eqnum,
                               const std::vector<double_complex> &psi   );

     /*! Returns the coefficients of the jacobian system associated with the
         Ginzburg--Landau equations. */
     void getJacobianRow( const int                         eqnum,
                          const std::vector<double_complex> &psi,
                          std::vector<int>                  &columnIndicesPsi,
                          std::vector<double_complex>       &valuesPsi,
                          std::vector<int>                  &columnIndicesPsiConj,
                          std::vector<double_complex>       &valuesPsiConj );

     /*! Get sparsity pattern of the jacobian system. */
     void getJacobianRowSparsity( int              eqnum,
                                  std::vector<int> &columnIndicesPsi,
                                  std::vector<int> &columnIndicesPsiConj );

     /*! Calcuate the grid approximation of the Gibbs free energy
       \f[
       \mathcal{G} = \int\limits_{\Omega} |\psi|^4 \,\mathrm{d}\omega
       \f]
       of a given state \f$\psi\f$. */
     double freeEnergy( const std::vector<double_complex> &psi );

  private:
      double h; //! mesh width
      PsiGrid::PsiGrid psiGrid;
      AGrid::AGrid     aGrid;

      void getEquationType( const int,
                            equationType&,
                            int* );

      enum filltype { VALUES, SPARSITY };

      /*! Calculated the coefficients of the jacobian system associated with the
          Ginzburg--Landau equations. */
      void computeJacobianRow( const filltype                    ft,
                               const int                         eqnum,
                               const std::vector<double_complex> &psi,
                               std::vector<int>                  &columnIndicesPsi,
                               std::vector<double_complex>       &valuesPsi,
                               std::vector<int>                  &columnIndicesPsiConj,
                               std::vector<double_complex>       &valuesPsiConj         );

};