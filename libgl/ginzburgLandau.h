/********************************************//**
 * The Ginzburg--Landau equations.
 ***********************************************/
#include <complex>

#include "psiGrid.h"
#include "aGrid.h"

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
     std::complex<double> computeGl( int eqnum,
                                     std::vector<std::complex<double> > psi );

     /*! Returns the coefficients of the jacobian system associated with the
         Ginzburg--Landau equations. */
     void getJacobianRow( int eqnum,
                          std::vector<std::complex<double> >  psi,
                          std::vector<int>&                   columnIndicesPsi,
                          std::vector<std::complex<double> >& valuesPsi,
                          std::vector<int>&                   columnIndicesPsiConj,
                          std::vector<std::complex<double> >& valuesPsiConj );

     /*! Get sparsity pattern of the jacobian system. */
     void getJacobianRowSparsity( int eqnum,
                                  std::vector<int>& columnIndicesPsi,
                                  std::vector<int>& columnIndicesPsiConj );

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
      void computeJacobianRow( filltype ft,
                               int eqnum,
                               std::vector<std::complex<double> > psi,
                               std::vector<int>& columnIndicesPsi,
                               std::vector<std::complex<double> >& valuesPsi,
                               std::vector<int>& columnIndicesPsiConj,
                               std::vector<std::complex<double> >& valuesPsiConj );

};