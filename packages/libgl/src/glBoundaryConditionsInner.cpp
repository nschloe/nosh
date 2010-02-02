#include "glBoundaryConditionsInner.h"

#include <Teuchos_Array.hpp>

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
GlBoundaryConditionsInner::getGlEntry ( const int                       eqIndex,
                                        const ComplexVector           & psi,
                                        const double                    chi,
                                        const GridUniformVirtual      & grid,
                                        const MagneticVectorPotential & A
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
      // -------------------------------------------------------------------
      psiKRight = psiView[ grid.getKRight(k) ];
      psiKAbove = psiView[ grid.getKAbove(k) ];

      xRight = grid.getXRight( k );
      xAbove = grid.getXAbove( k );

      ARight = A.getAx( *xRight );
      AAbove = A.getAy( *xAbove );

      res = ( - psiK      * 2.0
              + psiKRight * exp(-I*ARight*h)
              + psiKAbove * exp(-I*AAbove*h) ) * I/(sqrt(2)*h);
      res *= exp( I*chi );
      // -------------------------------------------------------------------
    break;

  case BOTTOMRIGHTCONVEX:
    // -----------------------------------------------------------------------
    psiKLeft  = psiView[ grid.getKLeft ( k ) ];
    psiKAbove = psiView[ grid.getKAbove ( k ) ];

    xLeft  = grid.getXLeft ( k );
    xAbove = grid.getXAbove( k );

    ALeft  = A.getAx( *xLeft );
    AAbove = A.getAy( *xAbove );

    res = ( - psiK      * 2.0
            + psiKLeft  * exp ( I*ALeft *h )
            + psiKAbove * exp ( -I*AAbove*h ) ) * I/ ( sqrt ( 2 ) *h );
    // -----------------------------------------------------------------------
    break;

  case TOPRIGHTCONVEX:
    // -----------------------------------------------------------------------
    psiKLeft  = psiView[ grid.getKLeft ( k ) ];
    psiKBelow = psiView[ grid.getKBelow ( k ) ];

    xLeft  = grid.getXLeft ( k );
    xBelow = grid.getXBelow( k );

    ALeft  = A.getAx( *xLeft );
    ABelow = A.getAy( *xBelow );

    res = ( - psiK      * 2.0
            + psiKLeft  * exp ( I*ALeft *h )
            + psiKBelow * exp ( I*ABelow*h ) ) * I/ ( sqrt ( 2 ) *h );
    res *= exp( I*chi );
    // -----------------------------------------------------------------------

    break;

  case TOPLEFTCONVEX:
    // -----------------------------------------------------------------------
    psiKRight = psiView[ grid.getKRight ( k ) ];
    psiKBelow = psiView[ grid.getKBelow ( k ) ];

    xRight = grid.getXRight( k );
    xBelow = grid.getXBelow( k );

    ARight = A.getAx( *xRight );
    ABelow = A.getAy( *xBelow );

    res = ( - psiK      * 2.0
            + psiKRight * exp ( -I*ARight*h )
            + psiKBelow * exp ( I*ABelow*h ) ) * I/ ( sqrt ( 2 ) *h );
    res *= exp( I*chi );
    // -----------------------------------------------------------------------
    break;

  case BOTTOM:
    // -------------------------------------------------------------------
    // normal derivative
    psiKAbove = psiView[ grid.getKAbove ( k ) ];

    xAbove = grid.getXAbove( k );
    AAbove = A.getAy( *xAbove );

    res = ( - psiK
            + psiKAbove * exp ( -I*AAbove*h ) ) * I/h;
    res *= exp( I*chi );
    // -------------------------------------------------------------------
    break;

  case RIGHT:
    // -------------------------------------------------------------------
    // normal derivative
    psiKLeft = psiView[ grid.getKLeft ( k ) ];

    xLeft  = grid.getXLeft ( k );
    ALeft  = A.getAx( *xLeft );

    res = ( - psiK
            + psiKLeft * exp ( I*ALeft*h ) ) * I/h;
    res *= exp( I*chi );
    // -------------------------------------------------------------------
    break;

  case TOP:
    // -------------------------------------------------------------------
    // normal derivative
    psiKBelow = psiView[ grid.getKBelow ( k ) ];

    xBelow = grid.getXBelow( k );
    ABelow = A.getAy( *xBelow );

    res = ( - psiK
            + psiKBelow * exp ( I*ABelow*h ) ) * I/h;
    res *= exp( I*chi );
    // -------------------------------------------------------------------
    break;

  case LEFT:
    // -------------------------------------------------------------------
    // normal derivative
    psiKRight = psiView[ grid.getKRight ( k ) ];

    xRight = grid.getXRight( k );
    ARight = A.getAx( *xRight );

    res = ( - psiK
            + psiKRight * exp ( -I*ARight*h ) ) * I/h;
    res *= exp( I*chi );
    // -------------------------------------------------------------------
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
GlBoundaryConditionsInner::getGlJacobianRow ( const int                           eqIndex,
                                              const Teuchos::RCP<ComplexVector> & psi,
                                              const double                        chi,
                                              const GridUniformVirtual          & grid,
                                              const MagneticVectorPotential     & A,
                                              const bool                          fillValues,
                                              std::vector<int>                  & columnIndicesPsi,
                                              std::vector<double_complex>       & valuesPsi,
                                              std::vector<int>                  & columnIndicesPsiConj,
                                              std::vector<double_complex>       & valuesPsiConj
                                            ) const
{
  int k, kLeft, kRight, kBelow, kAbove;
  int numEntriesPsi, numEntriesPsiConj;
  double ARight, ALeft, AAbove, ABelow;
  double h = grid.getUniformH();
  Teuchos::RCP<Teuchos::Array<double> > xRight = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xLeft  = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xAbove = Teuchos::rcp( new Teuchos::Array<double>(2) );
  Teuchos::RCP<Teuchos::Array<double> > xBelow = Teuchos::rcp( new Teuchos::Array<double>(2) );

  // associate the equation index with the grid point k
  k = eqIndex;

  equationType eqType;
  eqType = getEquationType ( k, grid );



  switch ( eqType )
    {
    case BOTTOMLEFTCONVEX:
          // -------------------------------------------------------------------
          kRight = grid.getKRight( k );
          kAbove = grid.getKAbove( k );

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kRight;
          columnIndicesPsi[2] = kAbove;

          if (fillValues) {
              xRight = grid.getXRight( k );
              xAbove = grid.getXAbove( k );
              ARight = A.getAx( *xRight );
              AAbove = A.getAy( *xAbove );

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

    case BOTTOMRIGHTCONVEX:
      // -----------------------------------------------------------------------
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

    case TOPRIGHTCONVEX:
      // -----------------------------------------------------------------------
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

    case TOPLEFTCONVEX:
      // -----------------------------------------------------------------------
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
      kAbove = grid.getKAbove ( k );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kAbove;

      if ( fillValues )
        {
          xAbove = grid.getXAbove( k );
          AAbove = A.getAy( *xAbove );

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
      kLeft = grid.getKLeft ( k );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;

      if ( fillValues )
        {
          xLeft = grid.getXLeft ( k );
          ALeft = A.getAx( *xLeft );

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
      kBelow = grid.getKBelow ( k );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;

      if ( fillValues )
        {
          xBelow = grid.getXBelow( k );
          ABelow = A.getAy( *xBelow );

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
      kRight = grid.getKRight ( k );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;

      if ( fillValues )
        {
          xRight = grid.getXRight( k );
          ARight = A.getAx( *xRight );

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
        TEST_FOR_EXCEPTION( psi.is_null(),
    		                std::logic_error,
    		                "Illegal equationType" << eqType  );
    }

}
// =============================================================================
