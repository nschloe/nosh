#include <complex>

#include "psiGrid.h"
#include "aGrid.h"

class GinzburgLandau
{
  public:
     // constructor
     GinzburgLandau( int nx );

     // destructor
     ~GinzburgLandau();

     // evaluate the 
     std::complex<double> boundaryConditions( int k,
                                              double* psiReal,
                                              double* psiImag,
                                              PsiGrid::PsiGrid,
                                              AGrid::AGrid );

  private:
      int Nx;
      int d;
      double h;
};