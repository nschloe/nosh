#include "glBoundaryConditionsOuter.h"

#include <Teuchos_Array.hpp>

// complex unit
const double_complex I ( 0,1 );

// =============================================================================
GlBoundaryConditionsOuter::GlBoundaryConditionsOuter()
{
}
// =============================================================================
GlBoundaryConditionsOuter::~GlBoundaryConditionsOuter()
{
}
// =============================================================================
double_complex
GlBoundaryConditionsOuter::getGlEntry ( const int                                eqIndex,
                                        const Tpetra::Vector<double_complex,int> &psi,
                                        const GridUniformVirtual                 &grid,
                                        const MagneticVectorPotential            &A
                                      ) const
{
  double_complex res;
  double_complex psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
  double ARight, ALeft, AAbove, ABelow;
  double h = grid.getUniformH();
  Teuchos::RCP<Teuchos::Array<double> > xRight = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xLeft  = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xAbove = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xBelow = Teuchos::rcp( new Teuchos::Array<double>(2) );

  // associate the equation index with the grid point k
  int k = eqIndex;

  equationType eqType;
  eqType = getEquationType ( k, grid );

  // Get a view of the whole vector.
  // Remember: This only works with one core.
  Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

  psiK = psiView[k];

  switch ( eqType )
    {
    case BOTTOMLEFTCONVEX:
      // -------------------------------------------------------------------------
      // interior equation, then outward derivative substituted
      psiKRight = psiView[ grid.getKRight ( k ) ];
      psiKAbove = psiView[ grid.getKAbove ( k ) ];

      xRight = grid.getXRight( k );
      xAbove = grid.getXAbove( k );

      ARight = A.getAx( *xRight );
      AAbove = A.getAy( *xAbove );

      res = ( - psiK      * 2.0
              + psiKRight * exp ( -I*ARight*h )
              + psiKAbove * exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // ------------------------------------------------------------------------
      break;

    case BOTTOMRIGHTCONVEX:
      // ---------------------------------------------------------------------------
      psiKLeft  = psiView[ grid.getKLeft ( k ) ];
      psiKAbove = psiView[ grid.getKAbove ( k ) ];

      xLeft  = grid.getXLeft( k );
      xAbove = grid.getXAbove( k );

      ALeft  = A.getAx( *xLeft );
      AAbove = A.getAy( *xAbove );

      res = ( psiK* ( -2.0 )
              + psiKLeft * exp ( I*ALeft *h )
              + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // ---------------------------------------------------------------------------
      break;

    case TOPRIGHTCONVEX:
      // ---------------------------------------------------------------------------
      psiKLeft  = psiView[ grid.getKLeft ( k ) ];
      psiKBelow = psiView[ grid.getKBelow ( k ) ];

      xLeft  = grid.getXLeft ( k );
      xBelow = grid.getXBelow( k );

      ALeft  = A.getAx( *xLeft );
      ABelow = A.getAy( *xBelow );

      res = ( psiK* ( -2.0 )
              + psiKLeft * exp ( I*ALeft *h )
              + psiKBelow* exp ( I*ABelow*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // ---------------------------------------------------------------------------

      break;

    case TOPLEFTCONVEX:
      // ---------------------------------------------------------------------------
      psiKRight = psiView[ grid.getKRight ( k ) ];
      psiKBelow = psiView[ grid.getKBelow ( k ) ];

      xRight = grid.getXRight( k );
      xBelow = grid.getXBelow( k );

      ARight = A.getAx( *xRight );
      ABelow = A.getAy( *xBelow );

      res = ( psiK* ( -2.0 )
              + psiKRight* exp ( -I*ARight*h )
              + psiKBelow* exp (  I*ABelow*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // ---------------------------------------------------------------------------
      break;

    case BOTTOM:
      // -----------------------------------------------------------------------
      psiKLeft  = psiView[ grid.getKLeft ( k ) ];
      psiKRight = psiView[ grid.getKRight ( k ) ];
      psiKAbove = psiView[ grid.getKAbove ( k ) ];

      xLeft  = grid.getXLeft( k );
      xRight = grid.getXRight( k );
      xAbove = grid.getXAbove( k );

      ALeft  = A.getAx( *xLeft );
      ARight = A.getAx( *xRight );
      AAbove = A.getAy( *xAbove );

      res = ( psiK* ( -3.0 )
              + psiKLeft*  exp ( I*ALeft *h ) + psiKRight* exp ( -I*ARight*h )
              + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1.0-norm ( psiK ) );
      // -----------------------------------------------------------------------
      break;

    case RIGHT:
      // -----------------------------------------------------------------------
      psiKLeft  = psiView[ grid.getKLeft ( k ) ];
      psiKBelow = psiView[ grid.getKBelow ( k ) ];
      psiKAbove = psiView[ grid.getKAbove ( k ) ];

      xLeft  = grid.getXLeft( k );
      xBelow = grid.getXBelow( k );
      xAbove = grid.getXAbove( k );

      ALeft  = A.getAx( *xLeft );
      ABelow = A.getAy( *xBelow );
      AAbove = A.getAy( *xAbove );

      res = ( psiK* ( -3.0 )
              + psiKLeft*  exp ( I*ALeft *h )
              + psiKBelow* exp ( I*ABelow*h ) + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // -----------------------------------------------------------------------
      break;

    case TOP:
      // -----------------------------------------------------------------------
      psiKLeft  = psiView[ grid.getKLeft ( k ) ];
      psiKRight = psiView[ grid.getKRight ( k ) ];
      psiKBelow = psiView[ grid.getKBelow ( k ) ];

      xLeft  = grid.getXLeft( k );
      xRight = grid.getXRight( k );
      xBelow = grid.getXBelow( k );

      ALeft  = A.getAx( *xLeft );
      ARight = A.getAx( *xRight );
      ABelow = A.getAy( *xBelow );

      res = ( psiK* ( -3.0 )
              + psiKLeft*  exp ( I*ALeft *h ) + psiKRight* exp ( -I*ARight*h )
              + psiKBelow* exp ( I*ABelow*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // -----------------------------------------------------------------------
      break;

    case LEFT:
      // -----------------------------------------------------------------------
      psiKRight = psiView[ grid.getKRight ( k ) ];
      psiKBelow = psiView[ grid.getKBelow ( k ) ];
      psiKAbove = psiView[ grid.getKAbove ( k ) ];

      xRight = grid.getXRight( k );
      xBelow = grid.getXBelow( k );
      xAbove = grid.getXAbove( k );

      ARight = A.getAx( *xRight );
      ABelow = A.getAy( *xBelow );
      AAbove = A.getAy( *xAbove );

      res = ( psiK* ( -3.0 )
              + psiKRight* exp ( -I*ARight*h )
              + psiKBelow* exp ( I*ABelow*h ) + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // -----------------------------------------------------------------------
      break;

    default:
        TEST_FOR_EXCEPTION( true,
    		            std::logic_error,
    			    "Illegal equationType " << eqType );
    }

  // return the result
  return res;
}
// =============================================================================
void
GlBoundaryConditionsOuter::getGlJacobianRow ( const int                                               eqIndex,
                                              const Teuchos::RCP<Tpetra::Vector<double_complex,int> > &psi,
                                              const GridUniformVirtual                                &grid,
                                              const MagneticVectorPotential                           &A,
                                              const bool                                              fillValues,
                                              std::vector<int>                                        &columnIndicesPsi,
                                              std::vector<double_complex>                             &valuesPsi,
                                              std::vector<int>                                        &columnIndicesPsiConj,
                                              std::vector<double_complex>                             &valuesPsiConj
                                            ) const
{
  int k, kLeft, kRight, kBelow, kAbove;
  int numEntriesPsi, numEntriesPsiConj;
  double ARight, ALeft, AAbove, ABelow;
  double h = grid.getUniformH();
  Teuchos::RCP<Teuchos::Array<int> >    i      = Teuchos::rcp( new Teuchos::Array<int>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xRight = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xLeft  = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xAbove = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xBelow = Teuchos::rcp( new Teuchos::Array<double>(2) );

  // associate the equation index with the grid point k
  k = eqIndex;

  equationType eqType;
  eqType = getEquationType ( k, grid );


  // Get a view of the whole vector.
  // Remember: This only works with one core.
  Teuchos::ArrayRCP<const double_complex> psiView;

  if ( fillValues )
    {
      TEST_FOR_EXCEPTION( !psi.is_valid_ptr(),
  		                  std::logic_error,
  		                  "Values are supposed to be filled in, but psi is invalid." );
      TEST_FOR_EXCEPTION( psi.is_null(),
  		                  std::logic_error,
  		                  "Values are supposed to be filled in, but psi is NULL." );
      psiView = psi->get1dView();
    }

  switch ( eqType )
    {
    case BOTTOMLEFTCONVEX:
      // -----------------------------------------------------------------------
      kRight = grid.getKRight ( k );
      kAbove = grid.getKAbove ( k );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;
      columnIndicesPsi[2] = kAbove;

      if ( fillValues )
        {
          xRight = grid.getXRight( k );
          xAbove = grid.getXAbove( k );

          ARight = A.getAx( *xRight );
          AAbove = A.getAy( *xAbove );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = - 2.0 / ( h*h )
                         + ( 1 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( -I*ARight*h ) / ( h*h );
          valuesPsi[2] = exp ( -I*AAbove*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;
      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
      // -----------------------------------------------------------------------
      break;

    case BOTTOMRIGHTCONVEX:
      // ---------------------------------------------------------------------------
      kLeft  = grid.getKLeft ( k );
      kAbove = grid.getKAbove ( k );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kAbove;

      if ( fillValues )
        {
          xLeft  = grid.getXLeft ( k );
          xAbove = grid.getXAbove( k );

          ALeft  = A.getAx( *xLeft );
          AAbove = A.getAy( *xAbove );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -2.0 / ( h*h )
                         + ( 1 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( I*ALeft *h ) / ( h*h );
          valuesPsi[2] = exp ( -I*AAbove*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;
      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
      // ---------------------------------------------------------------------------
      break;

    case TOPRIGHTCONVEX:
      // ---------------------------------------------------------------------------
      kLeft  = grid.getKLeft ( k );
      kBelow = grid.getKBelow ( k );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kBelow;

      if ( fillValues )
        {
          xLeft  = grid.getXLeft ( k );
          xBelow = grid.getXBelow( k );

          ALeft  = A.getAx( *xLeft );
          ABelow = A.getAy( *xBelow );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -2.0 / ( h*h )
                         + ( 1 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( I*ALeft *h ) / ( h*h );
          valuesPsi[2] = exp ( I*ABelow*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;
      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
      // ---------------------------------------------------------------------------
      break;

    case TOPLEFTCONVEX:
      // ---------------------------------------------------------------------------
      kRight = grid.getKRight ( k );
      kBelow = grid.getKBelow ( k );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;
      columnIndicesPsi[2] = kBelow;

      if ( fillValues )
        {
          xRight = grid.getXRight( k );
          xBelow = grid.getXBelow( k );

          ARight = A.getAx( *xRight );
          ABelow = A.getAy( *xBelow );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -2.0 / ( h*h )
                         + ( 1 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( -I*ARight*h ) / ( h*h );
          valuesPsi[2] = exp ( I*ABelow*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;
      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
      // ---------------------------------------------------------------------------
      break;

    case BOTTOM:
      // -----------------------------------------------------------------------
      // normal derivative
      kLeft  = grid.getKLeft ( k );
      kRight = grid.getKRight ( k );
      kAbove = grid.getKAbove ( k );

      numEntriesPsi = 4;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kRight;
      columnIndicesPsi[3] = kAbove;

      if ( fillValues )
        {
          xLeft  = grid.getXLeft ( k );
          xRight = grid.getXRight( k );
          xAbove = grid.getXAbove( k );

          ALeft  = A.getAx( *xLeft );
          ARight = A.getAx( *xRight );
          AAbove = A.getAy( *xAbove );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = - 3.0            / ( h*h )
                         + ( 1.0 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( I*ALeft *h ) / ( h*h );
          valuesPsi[2] = exp ( -I*ARight*h ) / ( h*h );
          valuesPsi[3] = exp ( -I*AAbove*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;

      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
      // -----------------------------------------------------------------------
      break;

    case RIGHT:
      // -----------------------------------------------------------------------
      kBelow = grid.getKBelow ( k );
      kAbove = grid.getKAbove ( k );
      kLeft  = grid.getKLeft ( k );

      numEntriesPsi = 4;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;
      columnIndicesPsi[2] = kAbove;
      columnIndicesPsi[3] = kLeft;

      if ( fillValues )
        {
          xLeft  = grid.getXLeft ( k );
          xBelow = grid.getXBelow( k );
          xAbove = grid.getXAbove( k );

          ALeft  = A.getAx( *xLeft );
          ABelow = A.getAy( *xBelow );
          AAbove = A.getAy( *xAbove );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = - 3.0            / ( h*h )
                         + ( 1.0 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( I*ABelow*h ) / ( h*h );
          valuesPsi[2] = exp ( -I*AAbove*h ) / ( h*h );
          valuesPsi[3] = exp ( I*ALeft *h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;

      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
      // -----------------------------------------------------------------------
      break;

    case TOP:
      // -----------------------------------------------------------------------
      kBelow = grid.getKBelow ( k );
      kRight = grid.getKRight ( k );
      kLeft  = grid.getKLeft ( k );

      numEntriesPsi = 4;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;
      columnIndicesPsi[2] = kLeft;
      columnIndicesPsi[3] = kRight;

      if ( fillValues )
        {
          xLeft  = grid.getXLeft ( k );
          xRight = grid.getXRight( k );
          xBelow = grid.getXBelow( k );

          ALeft  = A.getAx( *xLeft );
          ARight = A.getAx( *xRight );
          ABelow = A.getAy( *xBelow );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = - 3.0             / ( h*h )
                         + ( 1.0 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( I*ABelow*h ) / ( h*h );
          valuesPsi[2] = exp ( I*ALeft *h ) / ( h*h );
          valuesPsi[3] = exp ( -I*ARight*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;

      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
      // -----------------------------------------------------------------------
      break;

    case LEFT:
      // -----------------------------------------------------------------------
      kBelow = grid.getKBelow ( k );
      kAbove = grid.getKAbove ( k );
      kRight = grid.getKRight ( k );

      numEntriesPsi = 4;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;
      columnIndicesPsi[2] = kAbove;
      columnIndicesPsi[3] = kRight;

      if ( fillValues )
        {
          xRight = grid.getXRight( k );
          xBelow = grid.getXBelow( k );
          xAbove = grid.getXAbove( k );

          ARight = A.getAx( *xRight );
          ABelow = A.getAy( *xBelow );
          AAbove = A.getAy( *xAbove );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = - 3.0            / ( h*h )
                         + ( 1.0 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( I*ABelow*h ) / ( h*h );
          valuesPsi[2] = exp ( -I*AAbove*h ) / ( h*h );
          valuesPsi[3] = exp ( -I*ARight*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;

      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
      // -----------------------------------------------------------------------
      break;

    default:
        TEST_FOR_EXCEPTION( psi.is_null(),
    		                std::logic_error,
    		                "Illegal equationType" << eqType  );
    }

}
// =============================================================================
