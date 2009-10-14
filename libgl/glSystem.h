/*! Jacobian system for the Ginzburg--Landau problem.
 *  This routine can be used as an interface to NOX.
 ******************************************************************************/
#ifndef GLSYSTEM_H
#define GLSYSTEM_H

#include "ginzburgLandau.h"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ENull.hpp>

#include <Epetra_CrsMatrix.h>

#include <Teuchos_ParameterList.hpp>

#include <NOX_Epetra_Interface_Required.H> // NOX base class
#include <NOX_Epetra_Interface_Jacobian.H> // NOX base class

#include <LOCA_Epetra_Interface_Required.H> // LOCA base class

#include <NOX_Abstract_PrePostOperator.H>

#include <LOCA_Parameter_Vector.H>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

class GlSystem:
// public NOX::Epetra::Interface::Required,
      public NOX::Epetra::Interface::Jacobian,
      public LOCA::Epetra::Interface::Required
  {
  public:


    GlSystem ( GinzburgLandau::GinzburgLandau                               &gl,
               const Teuchos::RCP<Epetra_Comm>                              eComm,
               const bool                                                   &reverse, // Actually, this has nothing to do with the linear system.
                                                                                      // Get out, pls!
               const Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi=Teuchos::ENull()  );

    //! Destructor
    ~GlSystem();

    //! Evaluate the Ginzburg--Landau functions at a given state defined
    //! by the input vector x.
    bool computeF ( const Epetra_Vector &x,
                    Epetra_Vector       &F,
                    const NOX::Epetra::Interface::Required::FillType fillFlag = Residual );

    //! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
    //! at a given state defined by the input vector x.
    bool computeJacobian ( const Epetra_Vector &x,
                           Epetra_Operator     &Jac );

    //! Dummy preconditioner function. So far does nothing but throwing
    //! an exception when called.
    bool computePreconditioner ( const Epetra_Vector     &x,
                                 Epetra_Operator         &Prec,
                                 Teuchos::ParameterList* precParams=0 );

    //! Returns the current state. Not necessarily a solution to the problem!
    //! @return Reference-counted pointer to the current state.
    Teuchos::RCP<Epetra_Vector> getSolution();

    //! Returns the current Jacobian.
    //! @return Reference-counted pointer to the Jacobian.
    Teuchos::RCP<Epetra_CrsMatrix> getJacobian();

    //! Set the problem parameters.
    void setParameters ( const LOCA::ParameterVector &p );

    //! Set directory to where all output gets written.
    //! @param directory Name of the directory.
    void setOutputDir ( const string &directory );

    //! Print the solution x along with the continuation parameter conParam
    //! to a file. This function is called internally by LOCA to print
    //! solutions of each continuation step.
    void printSolution ( const Epetra_Vector &x,
                         double              conParam );

    //! Explictly print the solution x along with the problem parameters
    //! to the file fileName.
    void solutionToFile ( const Epetra_Vector          &x,
                          Teuchos::ParameterList &problemParams,
                          const std::string            &fileName );

  private:

    int realIndex2complexIndex ( const int realIndex );

    void real2complex ( const Epetra_Vector           &x,
                        vector<std::complex<double> > &psi );

    void real2psi ( const Epetra_Vector                     &x,
                    Tpetra::MultiVector<double_complex,int> &psi );

    void psi2real ( const Tpetra::MultiVector<double_complex,int> &psi,
                    Epetra_Vector                                 &x    );

    void complex2real ( const vector<double_complex> &psi,
                        Teuchos::RCP<Epetra_Vector>  realvec );

    void makeRealMap ( const Teuchos::RCP<const Tpetra::Map<int> >  complexMap,
                       Teuchos::RCP<Epetra_Map>               &realMap   );

    enum jacCreator { ONLY_GRAPH, VALUES };
    bool createJacobian ( const jacCreator    jc,
                          const Epetra_Vector &x );

    int NumRealUnknowns;
    int NumMyElements;
    int NumComplexUnknowns;

    bool reverse;
    
    GinzburgLandau::GinzburgLandau          Gl;
    const Teuchos::RCP<const Epetra_Comm>   EComm;
    Teuchos::RCP<const Teuchos::Comm<int> > TComm;
    Teuchos::RCP<Epetra_Map>                RealMap;
    Teuchos::RCP<Epetra_Map>                EverywhereMap;
    Teuchos::RCP<const Tpetra::Map<int> >   ComplexMap;
    Epetra_Vector                           *rhs;
    Teuchos::RCP<Epetra_CrsGraph>           Graph;
    Teuchos::RCP<Epetra_CrsMatrix>          jacobian;
    Teuchos::RCP<Epetra_Vector>             initialSolution;
    std::string                             outputDir;  //!< directory to where the output is written
  };

#endif // GLSYSTEM_H