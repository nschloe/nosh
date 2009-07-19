#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "ginzburgLandau.h"

class GlSystem
{
  public:
     GlSystem( int nx,
               GinzburgLandau::GinzburgLandau gl,
               Epetra_Comm& comm );

     ~GlSystem();

     bool computeF( const Epetra_Vector& x,
                    Epetra_Vector& FVec );

//   protected:
//      // Attributes visible to descendents
  private:
      enum complexPart { REAL, IMAGINARY };

      int  realIndex2psiIndex ( int realIndex );
      bool real2psi( Epetra_Vector realvec,
                     std::complex<double>* psi );

      int NumGlobalElements;
      int NumMyElements;
      int NumComplexUnknowns;
      GinzburgLandau::GinzburgLandau Gl;
      int Nx;
      Epetra_Comm* Comm;
      Epetra_Vector* rhs;
      Epetra_Map *StandardMap, 
                 *EverywhereMap;
};