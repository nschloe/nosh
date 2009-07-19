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
     GinzburgLandau( int nx,
                     double edgelength,
                     double h0 );

     // destructor
     ~GinzburgLandau();

     // evaluate the 
     std::complex<double> computeGl( int eqnum,
                                     std::complex<double>* psi );

     void computeJacobianBlocks( int eqnum,
                                 std::complex<double>* psi,
                                 int* columnIndicesPsi, 
                                 int* columnIndicesPsiConj,
                                 std::complex<double>* valuesPsi,
                                 std::complex<double>* valuesPsiConj );

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
};