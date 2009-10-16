#include "glBoundaryConditionsInner.h"

#include <EpetraExt_Utils.h> // for toString
#include <Teuchos_Array.hpp>

#include "glException.h"

// complex unit
const double_complex I ( 0,1 );

// =============================================================================
GlBoundaryConditionsInner::GlBoundaryConditionsInner()
{
}
// =============================================================================
GlBoundaryConditionsInner::~GlBoundaryConditionsInner()
{
}
// =============================================================================
double_complex
GlBoundaryConditionsInner::getGlEntry ( const int                                     eqIndex,
                                        const Tpetra::MultiVector<double_complex,int> &psi,
                                        const StaggeredGrid::StaggeredGrid            &sGrid
                                      )
{
  double_complex res;
  double_complex psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
  double ARight, ALeft, AAbove, ABelow;
  double h = sGrid.getH();
  Teuchos::Array<int> i(2);
  
  equationType eqType;
  getEquationType ( eqIndex,
                    sGrid,
                    eqType,
                    i       );

  // Get a view of the whole vector.
  // Remember: This only works with one core.
  Teuchos::ArrayRCP<const double_complex> psiView = psi.getVector(0)->get1dView();

  int k = sGrid.i2k ( i );
  psiK = psiView[k];
  
  switch ( eqType )
  {
    case BOTTOMLEFT:
      // -------------------------------------------------------------------
      psiKRight = psiView[ sGrid.getKRight(i) ];
      psiKAbove = psiView[ sGrid.getKAbove(i) ];

      ARight = sGrid.getAxRight( i );
      AAbove = sGrid.getAyAbove( i );

      res = ( - psiK      * 2.0
              + psiKRight * exp(-I*ARight*h)
              + psiKAbove * exp(-I*AAbove*h) ) * I/(sqrt(2)*h);
      // -------------------------------------------------------------------
    break;

  case BOTTOMRIGHT:
    // -----------------------------------------------------------------------
    psiKLeft  = psiView[ sGrid.getKLeft ( i ) ];
    psiKAbove = psiView[ sGrid.getKAbove ( i ) ];

    ALeft  = sGrid.getAxLeft ( i );
    AAbove = sGrid.getAyAbove ( i );

    res = ( - psiK      * 2.0
            + psiKLeft  * exp ( I*ALeft *h )
            + psiKAbove * exp ( -I*AAbove*h ) ) * I/ ( sqrt ( 2 ) *h );
    // -----------------------------------------------------------------------
    break;

  case TOPRIGHT:
    // -----------------------------------------------------------------------
    psiKLeft  = psiView[ sGrid.getKLeft ( i ) ];
    psiKBelow = psiView[ sGrid.getKBelow ( i ) ];

    ALeft  = sGrid.getAxLeft ( i );
    ABelow = sGrid.getAyBelow ( i );

    res = ( - psiK      * 2.0
            + psiKLeft  * exp ( I*ALeft *h )
            + psiKBelow * exp ( I*ABelow*h ) ) * I/ ( sqrt ( 2 ) *h );
    // -----------------------------------------------------------------------

    break;

  case TOPLEFT:
    // -----------------------------------------------------------------------
    psiKRight = psiView[ sGrid.getKRight ( i ) ];
    psiKBelow = psiView[ sGrid.getKBelow ( i ) ];

    ARight = sGrid.getAxRight ( i );
    ABelow = sGrid.getAyBelow ( i );

    res = ( - psiK      * 2.0
            + psiKRight * exp ( -I*ARight*h )
            + psiKBelow * exp ( I*ABelow*h ) ) * I/ ( sqrt ( 2 ) *h );
    // -----------------------------------------------------------------------
    break;

  case BOTTOM:
    // -------------------------------------------------------------------
    // normal derivative
    psiKAbove = psiView[ sGrid.getKAbove ( i ) ];
    AAbove = sGrid.getAyAbove ( i );
    res = ( - psiK
            + psiKAbove * exp ( -I*AAbove*h ) ) * I/h;
    // -------------------------------------------------------------------
    break;

  case RIGHT:
    // -------------------------------------------------------------------
    // normal derivative
    psiKLeft = psiView[ sGrid.getKLeft ( i ) ];
    ALeft  = sGrid.getAxLeft ( i );
    res = ( - psiK
            + psiKLeft * exp ( I*ALeft*h ) ) * I/h;
    // -------------------------------------------------------------------
    break;

  case TOP:
    // -------------------------------------------------------------------
    // normal derivative
    psiKBelow = psiView[ sGrid.getKBelow ( i ) ];
    ABelow = sGrid.getAyBelow ( i );
    res = ( - psiK
            + psiKBelow * exp ( I*ABelow*h ) ) * I/h;
    // -------------------------------------------------------------------
    break;

  case LEFT:
    // -------------------------------------------------------------------
    // normal derivative
    psiKRight = psiView[ sGrid.getKRight ( i ) ];
    ARight = sGrid.getAxRight ( i );
    res = ( - psiK
            + psiKRight * exp ( -I*ARight*h ) ) * I/h;
    // -------------------------------------------------------------------
    break;

  default:
    throw glException ( "GlBoundaryConditionsInner::getGlEntry",
                        "Illegal equationType "
                        + EpetraExt::toString ( eqType ) + "." );
  }

  // return the result
  return res;
}
// =============================================================================
void 
GlBoundaryConditionsInner::getGlJacobianRow ( const int                                                    eqIndex,
                                              const Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi,
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
  Teuchos::Array<int> i(2);

  equationType eqType;
  getEquationType ( eqIndex,
                    sGrid,
                    eqType,
                    i       );

  // needed everywhere
  k = sGrid.i2k ( i );

  switch ( eqType )
    {
    case BOTTOMLEFT:
          // -------------------------------------------------------------------
          kRight = sGrid.getKRight( i );
          kAbove = sGrid.getKAbove( i );

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kRight;
          columnIndicesPsi[2] = kAbove;

          if (fillValues) {
              ARight = sGrid.getAxRight( i );
              AAbove = sGrid.getAyAbove( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -2.0             * I/(sqrt(2)*h);
              valuesPsi[1] = exp(-I*ARight*h) * I/(sqrt(2)*h);
              valuesPsi[2] = exp(-I*AAbove*h) * I/(sqrt(2)*h);
          }

          numEntriesPsiConj = 0;
          columnIndicesPsiConj.resize(numEntriesPsiConj);
          if (fillValues)
              valuesPsiConj.resize(numEntriesPsiConj);
          // -------------------------------------------------------------------
      break;

    case BOTTOMRIGHT:
      // -----------------------------------------------------------------------
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
          valuesPsi[0] = -2.0                * I/ ( sqrt ( 2 ) *h );
          valuesPsi[1] = exp (  I*ALeft *h ) * I/ ( sqrt ( 2 ) *h );
          valuesPsi[2] = exp ( -I*AAbove*h ) * I/ ( sqrt ( 2 ) *h );
        }

      numEntriesPsiConj = 0;
      columnIndicesPsiConj.resize(numEntriesPsiConj);
      if (fillValues)
          valuesPsiConj.resize(numEntriesPsiConj);
      // -----------------------------------------------------------------------
      break;

    case TOPRIGHT:
      // -----------------------------------------------------------------------
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
          valuesPsi[0] = -2.0             * I/ ( sqrt ( 2 ) *h );
          valuesPsi[1] = exp ( I*ALeft *h ) * I/ ( sqrt ( 2 ) *h );
          valuesPsi[2] = exp ( I*ABelow*h ) * I/ ( sqrt ( 2 ) *h );
        }

      numEntriesPsiConj = 0;
      columnIndicesPsiConj.resize(numEntriesPsiConj);
      if (fillValues)
          valuesPsiConj.resize(numEntriesPsiConj);
      // -----------------------------------------------------------------------
      break;

    case TOPLEFT:
      // -----------------------------------------------------------------------
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
          valuesPsi[0] = -2.0                * I/ ( sqrt ( 2 ) *h );
          valuesPsi[1] = exp ( -I*ARight*h ) * I/ ( sqrt ( 2 ) *h );
          valuesPsi[2] = exp (  I*ABelow*h ) * I/ ( sqrt ( 2 ) *h );
        }

      numEntriesPsiConj = 0;
      columnIndicesPsiConj.resize(numEntriesPsiConj);
      if (fillValues)
          valuesPsiConj.resize(numEntriesPsiConj);
      // -----------------------------------------------------------------------
      break;

    case BOTTOM:
      // -------------------------------------------------------------------
      // normal derivative
      kAbove = sGrid.getKAbove ( i );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kAbove;

      if ( fillValues )
        {
          AAbove = sGrid.getAyAbove ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -1.0                * I/h;
          valuesPsi[1] = exp ( -I*AAbove*h ) * I/h;
        }

      numEntriesPsiConj = 0;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      if ( fillValues )
        valuesPsiConj.resize ( numEntriesPsiConj );
      // -------------------------------------------------------------------
      break;

    case RIGHT:
      // -------------------------------------------------------------------
      // normal derivative
      kLeft = sGrid.getKLeft ( i );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;

      if ( fillValues )
        {
          ALeft = sGrid.getAxLeft ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -1.0            * I/h;
          valuesPsi[1] = exp ( I*ALeft*h ) * I/h;
        }

      numEntriesPsiConj = 0;
      columnIndicesPsiConj.resize(numEntriesPsiConj);
      if (fillValues)
          valuesPsiConj.resize(numEntriesPsiConj);
      // -------------------------------------------------------------------
      break;

    case TOP:
      // -------------------------------------------------------------------
      // normal derivative
      kBelow = sGrid.getKBelow ( i );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;

      if ( fillValues )
        {
          ABelow = sGrid.getAyBelow ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -1.0             * I/h;
          valuesPsi[1] = exp ( I*ABelow*h ) * I/h;
        }

      numEntriesPsiConj = 0;
      columnIndicesPsiConj.resize(numEntriesPsiConj);
      if (fillValues)
          valuesPsiConj.resize(numEntriesPsiConj);
      // -------------------------------------------------------------------
      break;

    case LEFT:
      // -------------------------------------------------------------------
      // normal derivative
      kRight = sGrid.getKRight ( i );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;

      if ( fillValues )
        {
          ARight = sGrid.getAxRight ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -1.0                * I/h;
          valuesPsi[1] = exp ( -I*ARight*h ) * I/h;
        }

      numEntriesPsiConj = 0;
      columnIndicesPsiConj.resize(numEntriesPsiConj);
      if (fillValues)
          valuesPsiConj.resize(numEntriesPsiConj);
      // -------------------------------------------------------------------
      break;

    default:
      throw glException ( "GlBoundaryConditionsInner::getJacobianRow",
                          "Illegal equationType"
                          + EpetraExt::toString ( eqType ) + "." );
    }

}
// =============================================================================