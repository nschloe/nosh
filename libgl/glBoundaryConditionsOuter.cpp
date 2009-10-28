#include "glBoundaryConditionsOuter.h"

#include <EpetraExt_Utils.h> // for toString
#include <Teuchos_Array.hpp>

#include "glException.h"

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
GlBoundaryConditionsOuter::getGlEntry ( const int                                     eqIndex,
                                        const Tpetra::Vector<double_complex,int> &psi,
                                        const StaggeredGrid::StaggeredGrid            &sGrid
                                      )
{
  double_complex res;
  double_complex psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
  double ARight, ALeft, AAbove, ABelow;
  double h = sGrid.getH();
  Teuchos::Array<int> i ( 2 );

  equationType eqType;
  getEquationType ( eqIndex,
                    sGrid,
                    eqType,
                    i );

  // Get a view of the whole vector.
  // Remember: This only works with one core.
  Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

  int k = sGrid.i2k ( i );
  psiK = psiView[k];

  switch ( eqType )
    {
    case BOTTOMLEFT:
      // -------------------------------------------------------------------------
      // interior equation, then outward derivative substituted
      psiKRight = psiView[ sGrid.getKRight ( i ) ];
      psiKAbove = psiView[ sGrid.getKAbove ( i ) ];

      ARight = sGrid.getAxRight ( i );
      AAbove = sGrid.getAyAbove ( i );

      res = ( - psiK      * 2.0
              + psiKRight * exp ( -I*ARight*h )
              + psiKAbove * exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // ------------------------------------------------------------------------
      break;

    case BOTTOMRIGHT:
      // ---------------------------------------------------------------------------
      psiKLeft  = psiView[ sGrid.getKLeft ( i ) ];
      psiKAbove = psiView[ sGrid.getKAbove ( i ) ];

      ALeft  = sGrid.getAxLeft ( i );
      AAbove = sGrid.getAyAbove ( i );

      res = ( psiK* ( -2.0 )
              + psiKLeft * exp ( I*ALeft *h )
              + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // ---------------------------------------------------------------------------
      break;

    case TOPRIGHT:
      // ---------------------------------------------------------------------------
      psiKLeft  = psiView[ sGrid.getKLeft ( i ) ];
      psiKBelow = psiView[ sGrid.getKBelow ( i ) ];

      ALeft  = sGrid.getAxLeft ( i );
      ABelow = sGrid.getAyBelow ( i );

      res = ( psiK* ( -2.0 )
              + psiKLeft * exp ( I*ALeft *h )
              + psiKBelow* exp ( I*ABelow*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // ---------------------------------------------------------------------------

      break;

    case TOPLEFT:
      // ---------------------------------------------------------------------------
      psiKRight = psiView[ sGrid.getKRight ( i ) ];
      psiKBelow = psiView[ sGrid.getKBelow ( i ) ];

      ARight = sGrid.getAxRight ( i );
      ABelow = sGrid.getAyBelow ( i );

      res = ( psiK* ( -2.0 )
              + psiKRight* exp ( -I*ARight*h )
              + psiKBelow* exp (  I*ABelow*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // ---------------------------------------------------------------------------
      break;

    case BOTTOM:
      // -----------------------------------------------------------------------
      psiKLeft  = psiView[ sGrid.getKLeft ( i ) ];
      psiKRight = psiView[ sGrid.getKRight ( i ) ];
      psiKAbove = psiView[ sGrid.getKAbove ( i ) ];

      ALeft  = sGrid.getAxLeft ( i );
      ARight = sGrid.getAxRight ( i );
      AAbove = sGrid.getAyAbove ( i );
      res = ( psiK* ( -3.0 )
              + psiKLeft*  exp ( I*ALeft *h ) + psiKRight* exp ( -I*ARight*h )
              + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1.0-norm ( psiK ) );
      // -----------------------------------------------------------------------
      break;

    case RIGHT:
      // -----------------------------------------------------------------------
      psiKLeft  = psiView[ sGrid.getKLeft ( i ) ];
      psiKBelow = psiView[ sGrid.getKBelow ( i ) ];
      psiKAbove = psiView[ sGrid.getKAbove ( i ) ];

      ALeft  = sGrid.getAxLeft ( i );
      ABelow = sGrid.getAyBelow ( i );
      AAbove = sGrid.getAyAbove ( i );
      res = ( psiK* ( -3.0 )
              + psiKLeft*  exp ( I*ALeft *h )
              + psiKBelow* exp ( I*ABelow*h ) + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // -----------------------------------------------------------------------
      break;

    case TOP:
      // -----------------------------------------------------------------------
      psiKLeft  = psiView[ sGrid.getKLeft ( i ) ];
      psiKRight = psiView[ sGrid.getKRight ( i ) ];
      psiKBelow = psiView[ sGrid.getKBelow ( i ) ];

      ALeft  = sGrid.getAxLeft ( i );
      ARight = sGrid.getAxRight ( i );
      ABelow = sGrid.getAyBelow ( i );
      res = ( psiK* ( -3.0 )
              + psiKLeft*  exp ( I*ALeft *h ) + psiKRight* exp ( -I*ARight*h )
              + psiKBelow* exp ( I*ABelow*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // -----------------------------------------------------------------------
      break;

    case LEFT:
      // -----------------------------------------------------------------------
      psiKRight = psiView[ sGrid.getKRight ( i ) ];
      psiKBelow = psiView[ sGrid.getKBelow ( i ) ];
      psiKAbove = psiView[ sGrid.getKAbove ( i ) ];

      ARight = sGrid.getAxRight ( i );
      ABelow = sGrid.getAyBelow ( i );
      AAbove = sGrid.getAyAbove ( i );
      res = ( psiK* ( -3.0 )
              + psiKRight* exp ( -I*ARight*h )
              + psiKBelow* exp ( I*ABelow*h ) + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // -----------------------------------------------------------------------
      break;

    default:
      throw glException ( "GlBoundaryConditionsOuter::getGlEntry",
                          "Illegal equationType "
                          + EpetraExt::toString ( eqType ) + "." );
    }

  // return the result
  return res;
}
// =============================================================================
void
GlBoundaryConditionsOuter::getGlJacobianRow ( const int                                                    eqIndex,
    const Teuchos::RCP<Tpetra::Vector<double_complex,int> > psi,
    const StaggeredGrid::StaggeredGrid                           &sGrid,
    const bool                                                   fillValues,
    std::vector<int>                                             &columnIndicesPsi,
    std::vector<double_complex>                                  &valuesPsi,
    std::vector<int>                                             &columnIndicesPsiConj,
    std::vector<double_complex>                                  &valuesPsiConj
                                            )
{
  int k, kLeft, kRight, kBelow, kAbove;
  int numEntriesPsi, numEntriesPsiConj;
  double ARight, ALeft, AAbove, ABelow;
  double h = sGrid.getH();
  Teuchos::Array<int> i ( 2 );

  equationType eqType;
  getEquationType ( eqIndex,
                    sGrid,
                    eqType,
                    i );

  // needed everywhere
  k = sGrid.i2k ( i );

  // Get a view of the whole vector.
  // Remember: This only works with one core.
  Teuchos::ArrayRCP<const double_complex> psiView;

  if ( fillValues )
    {
      if ( !psi.is_valid_ptr() )
        {
          std::string message ( "Values are supposed to be filled in, but psi is invalid." );
          throw glException ( "GlBoundaryConditionsOuter::getGlJacobianRow",
                              message );
        }
      if ( psi.is_null() )
        {
          std::string message ( "Values are supposed to be filled in, but psi is NULL." );
          throw glException ( "GlBoundaryConditionsOuter::getGlJacobianRow",
                              message );
        }
      psiView = psi->get1dView();
    }

  switch ( eqType )
    {
    case BOTTOMLEFT:
      // -----------------------------------------------------------------------
      kRight = sGrid.getKRight ( i );
      kAbove = sGrid.getKAbove ( i );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;
      columnIndicesPsi[2] = kAbove;

      if ( fillValues )
        {
          ARight = sGrid.getAxRight ( i );
          AAbove = sGrid.getAyAbove ( i );

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

    case BOTTOMRIGHT:
      // ---------------------------------------------------------------------------
      kLeft  = sGrid.getKLeft ( i );
      kAbove = sGrid.getKAbove ( i );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kAbove;

      if ( fillValues )
        {
          ALeft    = sGrid.getAxLeft ( i );
          AAbove   = sGrid.getAyAbove ( i );
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

    case TOPRIGHT:
      // ---------------------------------------------------------------------------
      kLeft  = sGrid.getKLeft ( i );
      kBelow = sGrid.getKBelow ( i );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kBelow;

      if ( fillValues )
        {
          ALeft    = sGrid.getAxLeft ( i );
          ABelow   = sGrid.getAyBelow ( i );
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

    case TOPLEFT:
      // ---------------------------------------------------------------------------
      kRight = sGrid.getKRight ( i );
      kBelow = sGrid.getKBelow ( i );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;
      columnIndicesPsi[2] = kBelow;

      if ( fillValues )
        {
          ARight    = sGrid.getAxRight ( i );
          ABelow    = sGrid.getAyBelow ( i );
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
      kLeft  = sGrid.getKLeft ( i );
      kRight = sGrid.getKRight ( i );
      kAbove = sGrid.getKAbove ( i );

      numEntriesPsi = 4;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kRight;
      columnIndicesPsi[3] = kAbove;

      if ( fillValues )
        {
          ALeft  = sGrid.getAxLeft ( i );
          ARight = sGrid.getAxRight ( i );
          AAbove = sGrid.getAyAbove ( i );
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
      kBelow = sGrid.getKBelow ( i );
      kAbove = sGrid.getKAbove ( i );
      kLeft  = sGrid.getKLeft ( i );

      numEntriesPsi = 4;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;
      columnIndicesPsi[2] = kAbove;
      columnIndicesPsi[3] = kLeft;

      if ( fillValues )
        {
          ABelow = sGrid.getAyBelow ( i );
          AAbove = sGrid.getAyAbove ( i );
          ALeft  = sGrid.getAxLeft ( i );
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
      kBelow = sGrid.getKBelow ( i );
      kRight = sGrid.getKRight ( i );
      kLeft  = sGrid.getKLeft ( i );

      numEntriesPsi = 4;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;
      columnIndicesPsi[2] = kLeft;
      columnIndicesPsi[3] = kRight;

      if ( fillValues )
        {
          ABelow = sGrid.getAyBelow ( i );
          ALeft  = sGrid.getAxLeft ( i );
          ARight = sGrid.getAxRight ( i );
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
      kBelow = sGrid.getKBelow ( i );
      kAbove = sGrid.getKAbove ( i );
      kRight = sGrid.getKRight ( i );

      numEntriesPsi = 4;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;
      columnIndicesPsi[2] = kAbove;
      columnIndicesPsi[3] = kRight;

      if ( fillValues )
        {
          ABelow = sGrid.getAyBelow ( i );
          AAbove = sGrid.getAyAbove ( i );
          ARight = sGrid.getAxRight ( i );
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
      throw glException ( "GlBoundaryConditionsOuter::getJacobianRow",
                          "Illegal equationType"
                          + EpetraExt::toString ( eqType ) + "." );
    }

}
// =============================================================================
