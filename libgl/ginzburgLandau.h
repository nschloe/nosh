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
                                     std::complex<double>* psi );

     // get coefficients of the jacobian system
     void getJacobianRow( int eqnum,
                          std::complex<double>* psi,
                          int& numEntriesPsi,
                          int* columnIndicesPsi,
                          std::complex<double>* valuesPsi,
                          int& numEntriesPsiConj,
                          int* columnIndicesPsiConj,
                          std::complex<double>* valuesPsiConj );

     // get sparsity pattern of the jacobian system
     void getJacobianRowSparsity( int eqnum,
                                  int& numEntriesPsi,
                                  int* columnIndicesPsi,
                                  int& numEntriesPsiConj,
                                  int* columnIndicesPsiConj );

  private:
      int Nx;
      int d;
      double h;
      double Edgelength;
      double H0;
      PsiGrid::PsiGrid psiGrid;
      AGrid::AGrid     aGrid;

      void getEquationType( int,
                            equationType&,
                            int* );

      enum filltype { VALUES, SPARSITY };

      void computeJacobianRow( filltype ft,
                               int eqnum,
                               std::complex<double>* psi,
                               int& numEntriesPsi,
                               int* columnIndicesPsi,
                               std::complex<double>* valuesPsi,
                               int& numEntriesPsiConj,
                               int* columnIndicesPsiConj,
                               std::complex<double>* valuesPsiConj );
};