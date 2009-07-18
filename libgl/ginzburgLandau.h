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
     std::complex<double> boundaryConditions( int eqnum,
                                              double* psiReal,
                                              double* psiImag,
                                              PsiGrid::PsiGrid,
                                              AGrid::AGrid );

     std::complex<double> interiorEquations( int eqnum,
                                             std::complex<double>* psi,
                                             PsiGrid::PsiGrid psiGrid,
                                             AGrid::AGrid     aGrid );

  private:
      int Nx;
      int d;
      double h;
};