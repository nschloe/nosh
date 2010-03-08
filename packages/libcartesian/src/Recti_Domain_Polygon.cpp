/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "Recti_Domain_Polygon.h"

// ============================================================================
Recti::Domain::Polygon::
Polygon ( const Teuchos::Array<DoubleTuple> & polygonPoints ) :
        Recti::Domain::Abstract (),
        polygonPoints_ ( polygonPoints )
{
}
// ============================================================================
Recti::Domain::Polygon::
~Polygon()
{
}
// ============================================================================
Teuchos::Tuple<double,4>
Recti::Domain::Polygon::
getBoundingBox () const
{
    double xMin = polygonPoints_[0][0];
    double xMax = polygonPoints_[0][0];
    double yMin = polygonPoints_[0][1];
    double yMax = polygonPoints_[0][1];

    for ( int k=1; k<polygonPoints_.length(); k++ )
    {
        if ( polygonPoints_[k][0] < xMin )
            xMin = polygonPoints_[k][0];
        if ( polygonPoints_[k][0] > xMax )
            xMax = polygonPoints_[k][0];
        if ( polygonPoints_[k][1] < yMin )
            yMin = polygonPoints_[k][1];
        if ( polygonPoints_[k][1] > yMax )
            yMax = polygonPoints_[k][1];
    }

    return Teuchos::tuple ( xMin, yMin, xMax, yMax );
}
// ============================================================================
bool
Recti::Domain::Polygon::
isInDomain ( const DoubleTuple & x ) const
{
    // winding number algorithm
    // http://www.google.be/url?sa=t&source=web&ct=res&cd=1&ved=0CAkQFjAA&url=http%3A%2F%2Fwww.engr.colostate.edu%2F~dga%2Fdga%2Fpapers%2Fpoint_in_polygon.pdf&ei=ntxvS6eyOYLI-Qao2vjoAQ&usg=AFQjCNF_btLpRhTfkUQt34nXfsRylaF95g&sig2=Hh8I-SjpWC6queaj3QMeZw
    // TODO make w INT
    double w = 0;

    const double tol= 1.0e-15;

    Teuchos::RCP<DoubleTuple> currentNode;
    Teuchos::RCP<DoubleTuple> nextNode = Teuchos::rcp ( new DoubleTuple ( Teuchos::tuple ( polygonPoints_[0][0]-x[0], polygonPoints_[0][1]-x[1] ) ) );
    for ( int k=0; k<polygonPoints_.length(); k++ )
    {
        currentNode = nextNode;
        if ( k<polygonPoints_.length()-1 )
            nextNode = Teuchos::rcp ( new DoubleTuple ( Teuchos::tuple ( polygonPoints_[k+1][0]-x[0], polygonPoints_[k+1][1]-x[1] ) ) );
        else
            nextNode = Teuchos::rcp ( new DoubleTuple ( Teuchos::tuple ( polygonPoints_[0][0]-x[0], polygonPoints_[0][1]-x[1] ) ) );

        bool isK0Neg = ( *currentNode ) [1] < 0;
        bool isK1Neg = ( *nextNode ) [1] < 0;
        // check crossings of the positive x-axis
        if ( isK0Neg ^ isK1Neg ) // polygon edge crosses real axis
        {
            // calculate the intersection point
            double r = ( *currentNode ) [0] + ( *currentNode ) [1] * ( ( *nextNode ) [0]- ( *currentNode ) [0] ) / ( ( *currentNode ) [1]- ( *nextNode ) [1] );
            if ( r>0.0 )
            {
                if ( isK0Neg ) // crossing counter-clockwise
                    ++w;
                else  // crossing clockwise
                    --w;
            }
            else if ( fabs ( r ) <tol )
            { // point sits on the intersection of currentNode and nextNode
                return true;
            }
        }
        else if ( fabs ( ( *currentNode ) [1] ) < tol )
        {
            if ( ( *currentNode ) [0] > 0.0 )
            {
                if ( !isK1Neg )
                    w += 0.5;
                else
                    w-= 0.5;
            }
            else if ( fabs ( ( *currentNode ) [0] ) < tol     // point sits in currentNode
                      || fabs ( ( *nextNode ) [1] ) < tol ) // point sits horizontally between currentNode and nextNode
            {
                return true;
            }
        }
        else if ( fabs ( ( *nextNode ) [1] ) < tol && ( *nextNode ) [0] > 0.0 )
        {
            if ( isK0Neg )
                w += 0.5;
            else
                w-= 0.5;
        }
    }

    return w>0;
}
// ============================================================================
