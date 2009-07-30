/*! Jacobian system for the Ginzburg--Landau problem.
 *  This routine can be used as an interface to NOX.
 ******************************************************************************/
#include "ginzburgLandau.h"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>

#include <Epetra_CrsMatrix.h>

#include <Teuchos_ParameterList.hpp>

#include <NOX_Epetra_Interface_Required.H> // NOX base class
#include <NOX_Epetra_Interface_Jacobian.H> // NOX base class

#include <LOCA_Epetra_Interface_Required.H> // LOCA base class

#include <LOCA_Parameter_Vector.H>

class GlSystem: public NOX ::Epetra::Interface::Required,
                public NOX ::Epetra::Interface::Jacobian
//                 public LOCA::Epetra::Interface::Required
{
  public:

      //! Default constructor. 
      GlSystem( GinzburgLandau::GinzburgLandau &gl,
                Epetra_Comm                    &comm,
                std::vector<double_complex>    *psi = NULL );

      //! Destructor
      ~GlSystem();

      //! Evaluate the Ginzburg--Landau functions at a given state defined
      //! by the input vector x.
      bool computeF( const Epetra_Vector &x,
                     Epetra_Vector       &F,
                     const NOX::Epetra::Interface::Required::FillType fillFlag = Residual  );

      //! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
      //! at a given state defined by the input vector x.
      bool computeJacobian ( const Epetra_Vector &x,
                             Epetra_Operator     &Jac );

      //! Dummy preconditioner function. So far does nothing but throwing
      //! an exception when called.
      bool computePreconditioner( const Epetra_Vector     &x,
                                  Epetra_Operator         &Prec,
                                  Teuchos::ParameterList* precParams=0 );

      //! Returns the current state. Not necessarily a solution to the problem!
      Teuchos::RCP<Epetra_Vector> getSolution();

      //! Returns the current Jacobian.
      Teuchos::RCP<Epetra_CrsMatrix> getJacobian();

      void setParameters(const LOCA::ParameterVector &p);

      //! Set directory to where all output gets written.
      void setOutputDir( const string &directory );

      void printSolution( const Epetra_Vector &x,
                          double              conParam );

      void solutionToFile( const Epetra_Vector          &x,
                           const Teuchos::ParameterList &problemParams,
                           const std::string            &fileFormat,
                           const std::string            &fileName       );

  private:

      int realIndex2complexIndex ( const int realIndex );

      void real2complex( const Epetra_Vector           &x,
                         vector<std::complex<double> > &psi );

      void complex2real( const vector<double_complex> &psi,
                         Teuchos::RCP<Epetra_Vector>  realvec );

      bool initializeSoln();

      enum jacCreator { ONLY_GRAPH, VALUES };
      bool createJacobian( const jacCreator    jc,
                           const Epetra_Vector &x );

      int NumGlobalElements,
          NumMyElements,
          NumComplexUnknowns;

      GinzburgLandau::GinzburgLandau Gl;
      Epetra_Comm                    *Comm;
      Epetra_Map                     *StandardMap,
                                     *EverywhereMap;
      Epetra_Vector                  *rhs;
      Epetra_CrsGraph                *Graph;
      Teuchos::RCP<Epetra_CrsMatrix> jacobian;
      Teuchos::RCP<Epetra_Vector>    initialSolution;
      std::string                    outputDir;  //!< directory to where the output is written
};