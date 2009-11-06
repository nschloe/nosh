#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Comm.h>

#include <Tpetra_Vector.hpp>

#include <NOX.H>
#include <NOX_Epetra.H>

#include "Grid.h"
#include "glSystem.h"

// abbreviate the complex type name
typedef std::complex<double> double_complex;

class glNox
{

  public:

      // Class constructor 1
      glNox( const std::string fileName,
             const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
             const Teuchos::RCP<const Epetra_Comm>         &eComm );

      // Class constructor 2
      glNox( const int Nx,
             const double edgeLength,
             const double H0,
             const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
             const Teuchos::RCP<const Epetra_Comm>         &eComm );

      void
      setVerbose( bool verbose );

      void
      setSolverOptions( int maxNonlinearIterations );

      void
      createSolver();

      void
      createSolverGroup();

      void
      createConvergenceTests();

      void
      solve();

      double
      computeJacobianConditionNumber();

      void
      computeJacobianEigenvalues();

      void
      printSolutionToFile( std::string fileName );

      int
      checkConvergence();

  protected:

  private:
      const Teuchos::RCP<const Teuchos::Comm<int> > Comm_;
      const Teuchos::RCP<const Epetra_Comm>         eComm_;
      Teuchos::ParameterList problemParameters_;
      int MyPID_;
      Teuchos::RCP<GlSystem> glSystem_;
      Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr_;
      Teuchos::RCP<NOX::Solver::Generic> solver_;
      Teuchos::RCP<NOX::StatusTest::Combo> combo_;
      bool verbose_;
      Teuchos::RCP<NOX::Epetra::Group> grpPtr_;
      int maxNonlinearIterations_;

      void
      initializeGrid();

      void
      reOrder( Tpetra::Vector<double_complex> &psi,
               const Teuchos::RCP<Grid>       &grid );

      void
      setNonlinearSolverParameters( Teuchos::ParameterList & nlParams );

      void
      setPrintParameters( Teuchos::ParameterList & printParams );

      void
      setDirectionParameters( Teuchos::ParameterList & dirParams );

      void
      setSearchParameters( Teuchos::ParameterList & searchParams );

      void
      setNewtonParameters( Teuchos::ParameterList & newtonParams );

      void
      setLinearSolverParameters( Teuchos::ParameterList & lsParams );

};
