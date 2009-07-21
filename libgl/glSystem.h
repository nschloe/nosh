#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "ginzburgLandau.h"
#include "Teuchos_RCP.hpp"

#include "Epetra_CrsMatrix.h"

#include "Teuchos_ParameterList.hpp"

#include "NOX_Epetra_Interface_Required.H" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class
#include "NOX_Epetra_Interface_Preconditioner.H" // base class

class GlSystem: public NOX::Epetra::Interface::Required,
                public NOX::Epetra::Interface::Jacobian
{
  public:
     GlSystem( int nx,
               double h0,
               double edgelength,
               Epetra_Comm& comm );

     ~GlSystem();

     bool computeF( const Epetra_Vector& x,
                    Epetra_Vector& F,
                    const NOX::Epetra::Interface::Required::FillType fillFlag = Residual  );

     bool computeJacobian ( const Epetra_Vector &x,
                            Epetra_Operator &Jac    );

     bool computePreconditioner( const Epetra_Vector& x,
                                 Epetra_Operator& Prec,
                                 Teuchos::ParameterList* precParams=0 );

     Teuchos::RCP<Epetra_Vector>    getSolution();
     Teuchos::RCP<Epetra_CrsMatrix> getJacobian();

  private:
      int  realIndex2complexIndex ( int realIndex );

      void real2psi( Epetra_Vector realvec,
                     vector<std::complex<double> > psi );

      bool initializeSoln();

      bool createGraph();

      int NumGlobalElements;
      int NumMyElements;
      int NumComplexUnknowns;
      GinzburgLandau::GinzburgLandau Gl;
      int Nx;
      double H0;
      double Edgelength;
      Epetra_Comm* Comm;
      Epetra_Map *StandardMap, 
                 *EverywhereMap;
      Epetra_Vector* rhs;
      Epetra_CrsGraph* Graph;

      Teuchos::RCP<Epetra_CrsMatrix> jacobian;
      Teuchos::RCP<Epetra_Vector>    initialSolution;
};