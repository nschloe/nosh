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
        INTERIOR
      };

     // constructor
     GinzburgLandau( int nx );

     // destructor
     ~GinzburgLandau();

     // evaluate the 
     std::complex<double> computeGl( int eqnum,
                                     std::complex<double>* psi,
                                     PsiGrid::PsiGrid,
                                     AGrid::AGrid );

     void computeJacobianBlocks( int eqnum,
                                 std::complex<double>* psi,
                                 PsiGrid::PsiGrid psiGrid,
                                 AGrid::AGrid     aGrid,
                                 int* columnIndicesPsi, 
                                 int* columnIndicesPsiConj,
                                 std::complex<double>* valuesPsi,
                                 std::complex<double>* valuesPsiConj );

  private:
      int Nx;
      int d;
      double h;

      void getEquationType( int,
                            equationType&,
                            int* );
};