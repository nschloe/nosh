/*
 * GridVirtual.cpp
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schl\"omer
 */

#include "Recti_Grid_Abstract.h"

// ============================================================================
Recti::Grid::Abstract::
Abstract ( const Point        & h,
           const double         gridDomainArea,
           const unsigned int   numGridPoints,
           const double         scaling ) :
        h_ ( h ),
        scaling_ ( scaling ),
        gridDomainArea_ ( gridDomainArea ),
        numGridPoints_ ( numGridPoints )
{
}
// ============================================================================
Recti::Grid::Abstract::
Abstract() :
        h_ ( Teuchos::tuple(0.0,0.0,0.0) ),
        scaling_ ( 1.0 ),
        gridDomainArea_ ( 0.0 ),
        numGridPoints_ ( 0 )
{
}
// ============================================================================
Recti::Grid::Abstract::
~Abstract()
{
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Recti::Grid::Abstract::
getParameters() const
{
    Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp ( new LOCA::ParameterVector() );
  
    p->addParameter( "scaling", scaling_ );
    return p;
}
// ============================================================================
void
Recti::Grid::Abstract::
updateScaling ( const double scaling )
{    
    TEUCHOS_ASSERT_INEQUALITY ( scaling, >, 0.0 );

    double ratio = scaling/scaling_;

    for ( int k=0; k<h_.size(); k++ )
        h_[k] *= ratio;
    gridDomainArea_ *= ratio*ratio;  // rescale domain area
    scaling_         = scaling;
    
    return;
}
// ============================================================================
void
Recti::Grid::Abstract::
updateScaling ( const LOCA::ParameterVector & p )
{    
    if ( p.isParameter("scaling") )
    {
        double scaling = p.getValue ( "scaling" );
        updateScaling( scaling );
    }
    return;
}
// ============================================================================
Point
Recti::Grid::Abstract::
getH() const
{
    return h_;
}
// ============================================================================
unsigned int
Recti::Grid::Abstract::
getNumGridPoints() const
{
    return numGridPoints_;
}
// ============================================================================
double
Recti::Grid::Abstract::
getGridDomainArea() const
{
    return gridDomainArea_;
}
// =============================================================================
