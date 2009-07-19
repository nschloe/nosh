#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "ginzburgLandau.h"
#include "Teuchos_RCP.hpp"

#include "Epetra_CrsMatrix.h"

class GlSystem
{
  public:
     GlSystem( int nx,
               GinzburgLandau::GinzburgLandau gl,
               Epetra_Comm& comm );

     ~GlSystem();

     bool computeF( const Epetra_Vector& x,
                    Epetra_Vector& FVec );

     bool computeJacobian( const Epetra_Vector& x );

  private:
      enum complexPart { REAL, IMAGINARY };

      int  realIndex2psiIndex ( int realIndex );
      void real2psi( Epetra_Vector realvec,
                     std::complex<double>* psi );

      int NumGlobalElements;
      int NumMyElements;
      int NumComplexUnknowns;
      GinzburgLandau::GinzburgLandau Gl;
      int Nx;
      Epetra_Comm* Comm;
      Epetra_Map *StandardMap, 
                 *EverywhereMap;
      Epetra_Vector* rhs;
      Teuchos::RCP<Epetra_CrsMatrix> jacobian;
};