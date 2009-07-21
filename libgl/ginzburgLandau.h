#include <complex>

#include "psiGrid.h"
#include "aGrid.h"

class GinzburgLandau
{
  public:
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

     // constructor
     GinzburgLandau( int nx,
                     double edgelength,
                     double h0 );

     // destructor
     ~GinzburgLandau();

     // evaluate the 
     std::complex<double> computeGl( int eqnum,
                                     std::vector<std::complex<double> > psi );

     // get coefficients of the jacobian system
     void getJacobianRow( int eqnum,
                          std::vector<std::complex<double> > psi,
                          std::vector<int>& columnIndicesPsi,
                          std::vector<std::complex<double> >& valuesPsi,
                          std::vector<int>& columnIndicesPsiConj,
                          std::vector<std::complex<double> >& valuesPsiConj );

     // get sparsity pattern of the jacobian system
     void getJacobianRowSparsity( int eqnum,
                                  std::vector<int>& columnIndicesPsi,
                                  std::vector<int>& columnIndicesPsiConj );

  private:
      double h;
      PsiGrid::PsiGrid psiGrid;
      AGrid::AGrid     aGrid;

      void getEquationType( int,
                            equationType&,
                            int* );

      enum filltype { VALUES, SPARSITY };

      void computeJacobianRow( filltype ft,
                               int eqnum,
                               std::vector<std::complex<double> > psi,
                               std::vector<int>& columnIndicesPsi,
                               std::vector<std::complex<double> >& valuesPsi,
                               std::vector<int>& columnIndicesPsiConj,
                               std::vector<std::complex<double> >& valuesPsiConj );

};