// see <http://old.nabble.com/Undefined-reference-to-%27main%27-with-Boost-Test.-Why--td15986217.html>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GRNN operatorTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <Teuchos_DefaultComm.hpp>

#include "Ginla_FDM_State.h"
#include "Ginla_FDM_Operator_BCCentral.h"
#include "Ginla_MagneticVectorPotential_Centered.h"

#include "Recti_Domain_Circle.h"
#include "Recti_Domain_Square.h"

#include "Recti_Grid_Uniform.h"
#include "Recti_Grid_Reader.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// ===========================================================================
// <http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/tutorials/hello-the-testing-world.htmlhh>
BOOST_AUTO_TEST_CASE( operator_test )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( NULL, NULL );
#endif

    // Create a communicator for Tpetra objects
    const Teuchos::RCP<const Teuchos::Comm<int> > Comm =
            Teuchos::DefaultComm<int>::getComm();

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
            Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm>  eComm =
            Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

    // ------------------------------------------------------------------------
    // Create grid that contains all node types.
    // Let Omega be the unit circle.
    Teuchos::RCP<Recti::Domain::Abstract> Omega
         = Teuchos::rcp( new Recti::Domain::Circle(1.0) );

    double h = 5.0e-2;
    
    Teuchos::RCP<Recti::Grid::Uniform> grid
         = Teuchos::rcp( new Recti::Grid::Uniform( Omega, h ) );
    // ------------------------------------------------------------------------
    // create state on the grid
    Teuchos::RCP<Ginla::FDM::State> state
        = Teuchos::rcp( new Ginla::FDM::State( Comm, grid ) );
    state->getPsiNonConst()->putScalar( 1.0 );
    // ------------------------------------------------------------------------
    // create operator
    double h0 = 1.0;
    double edgeLength = 1.0;
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A
        = Teuchos::rcp( new Ginla::MagneticVectorPotential::Centered( h0 ) );
                                                                      
    Teuchos::RCP<Ginla::FDM::Operator::BCCentral> glOperator
        = Teuchos::rcp( new Ginla::FDM::Operator::BCCentral( grid,
                                                             A,
                                                             state->getPsi()->getMap(),
                                                             state->getPsi()->getMap() ) );
    // ------------------------------------------------------------------------
    // get dF/dh0 at state
    Teuchos::RCP<const Ginla::State::Virtual> dFState = glOperator->getDFDh0( state );
    // ------------------------------------------------------------------------
    // compute the finite difference
    double eps = 1.0e-05;
    LOCA::ParameterVector p;
    p.addParameter( "H0", h0 );
    
    p[0] = h0 + eps;
    glOperator->setParameters( p );
    Teuchos::RCP<Ginla::State::Virtual> statePlus = glOperator->getF( state );
    
    p[0] = h0 - eps;
    glOperator->setParameters( p );
    Teuchos::RCP<Ginla::State::Virtual> stateMinus = glOperator->getF( state );
    
    // store the finite difference in statePlus
    statePlus->update( -0.5/eps, *stateMinus, 0.5/eps );
    // ------------------------------------------------------------------------
    // check that the finite difference is somewhere near the analytic expression
    statePlus->update( -1.0, *dFState, 1.0 );
    
    BOOST_CHECK_SMALL( statePlus->normalizedScaledL2Norm(), 1.0e-05 );
    // ------------------------------------------------------------------------
    // clean up
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    // ------------------------------------------------------------------------

    return;
}
// ============================================================================
