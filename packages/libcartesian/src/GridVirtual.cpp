/*
 * GridVirtual.cpp
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schl\"omer
 */

#include "GridVirtual.h"

// ============================================================================
GridVirtual::GridVirtual ( const DoubleTuple  & h,
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
GridVirtual::GridVirtual() :
        h_ ( Teuchos::tuple<double> ( 0.0,0.0 ) ),
        scaling_ ( 1.0 ),
        gridDomainArea_ ( 0.0 ),
        numGridPoints_ ( 0 )
{
}
// ============================================================================
GridVirtual::~GridVirtual()
{
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
GridVirtual::getParameters() const
{
    Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp ( new LOCA::ParameterVector() );
  
    p->addParameter( "scaling", scaling_ );
    return p;
}
// ============================================================================
void
GridVirtual::updateScaling ( const double scaling )
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
GridVirtual::updateScaling ( const LOCA::ParameterVector & p)
{    
    TEST_FOR_EXCEPTION ( !p.isParameter ( "scaling" ),
                         std::logic_error,
                         "Label \"scaling\" not valid." );
    double scaling = p.getValue ( "scaling" );
    updateScaling(scaling);
    return;
}
// ============================================================================
DoubleTuple
GridVirtual::getH() const
{
    return h_;
}
// ============================================================================
unsigned int
GridVirtual::getNumGridPoints() const
{
    return numGridPoints_;
}
// ============================================================================
double
GridVirtual::getGridDomainArea() const
{
    return gridDomainArea_;
}
// =============================================================================
