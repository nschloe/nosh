/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

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

#include "Ginla_Helpers.h"

#include "NOX_Abstract_Vector.H"
#include "NOX_Epetra_Vector.H"

#include <Tpetra_Vector.hpp>

typedef std::complex<double> double_complex;

// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::Helpers::
teuchosParameterList2locaParameterVector( const Teuchos::ParameterList & p
                                        )
{
  Teuchos::RCP<LOCA::ParameterVector> pL =
      Teuchos::rcp( new LOCA::ParameterVector() );
        
  Teuchos::ParameterList::ConstIterator k;
  double * dummy;
  for ( k=p.begin(); k!=p.end(); ++k )
  {
    Teuchos::ParameterEntry e = p.entry(k);
    if ( e.isType<double>() )
      pL->addParameter( p.name(k),
                        e.getValue<double>(dummy) );
  }
  
  return pL;
}
// ============================================================================
Teuchos::RCP<Teuchos::ParameterList>
Ginla::Helpers::
locaParameterVector2teuchosParameterList( const LOCA::ParameterVector & pL )
{

    Teuchos::RCP<Teuchos::ParameterList> p =
      Teuchos::rcp( new Teuchos::ParameterList() );
      
    appendToTeuchosParameterList(  *p, pL );

    return p;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::Helpers::
mergeLocaParameterVectors( const LOCA::ParameterVector & p0,
                           const LOCA::ParameterVector & p1
                         )
{
  // intialize p with p0
  Teuchos::RCP<LOCA::ParameterVector> p = 
          Teuchos::rcp( new LOCA::ParameterVector( p0 ) );

  // add elements from p1
  for( int k=0; k<p1.length(); k++ )
  {
    double value = p1.getValue(k);
    std::string label = p1.getLabel(k);
    if( p->isParameter(label) )
    {
        // If the entry already exists, make sure the values
        // coincide.
        TEUCHOS_ASSERT_EQUALITY( p->getValue(label), value );
    }
    else
    {
        p->addParameter( label, value );
    }
  }
  
  return p;
}
// ============================================================================
void
Ginla::Helpers::
appendToTeuchosParameterList( Teuchos::ParameterList      & p,
                              const LOCA::ParameterVector & pL,
                              const std::string           & labelPrepend
                            )
{       
  Teuchos::ParameterList::ConstIterator k;
  for ( int k=0; k<pL.length(); k++ )
     p.set<double>( labelPrepend + pL.getLabel(k), pL[k] );  
  
  return;
}
// =============================================================================
// calculate the free energy of a state
double
Ginla::Helpers::
freeEnergy ( const ComplexVector        & psi,
             const Recti::Grid::General & grid )
{
    double localEnergy = 0.0;

    Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

    // sum up the energy on each processor
    for ( unsigned int k=0; k<psi.getLocalLength(); k++ )
    {
        int kGlobal = psi.getMap()->getGlobalElement ( k );
        double area = grid.cellArea ( kGlobal );
        localEnergy -= area * pow ( norm ( psiView[k] ),2 );
    }

    // reduce and scatter such that energy is available on
    // all cores
    int count = 1; // send *one* integer
    int numProcs =  psi.getMap()->getComm()->getSize();
    Teuchos::Array<double> sendBuff ( count ), recvBuff ( count );

    // fill send buffer
    sendBuff[0] = localEnergy;

    Teuchos::Array<int> recvCounts ( numProcs );
// fill recvCounts with {1,...,1}
    int numItemsPerProcess = 1;
    std::fill ( recvCounts.begin(), recvCounts.end(), numItemsPerProcess );

    Teuchos::reduceAllAndScatter ( * ( psi.getMap()->getComm() ),
                                   Teuchos::REDUCE_SUM,
                                   count,
                                   &sendBuff[0],
                                   &recvCounts[0],
                                   &recvBuff[0]
                                 );

    // normalize
    double maxEnergy = 1.0 * grid.getGridDomainArea();
    double globalEnergy = recvBuff[0] / maxEnergy;

    return globalEnergy;
}
// =============================================================================
double
Ginla::Helpers::
normalizedScaledL2Norm ( const ComplexVector        & psi,
                         const Recti::Grid::General & grid )
{
    double localSum = 0.0;

    Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

    // sum up on each processor
    for ( unsigned int k=0; k<psi.getLocalLength(); k++ )
    {
        int kGlobal = psi.getMap()->getGlobalElement ( k );
        double area = grid.cellArea ( kGlobal );
        localSum += area * norm ( psiView[k] );
    }

    // reduce and scatter such that energy is available on
    // all cores
    int count = 1; // send *one* integer
    Teuchos::Array<double> sendBuff ( count ), recvBuff ( count );

    // fill send buffer
    sendBuff[0] = localSum;

    int numProcs =  psi.getMap()->getComm()->getSize();
    Teuchos::Array<int> recvCounts ( numProcs );
// fill recvCounts with {1,...,1}
    int numItemsPerProcess = 1;
    std::fill ( recvCounts.begin(), recvCounts.end(), numItemsPerProcess );

    Teuchos::reduceAllAndScatter ( * ( psi.getMap()->getComm() ),
                                   Teuchos::REDUCE_SUM,
                                   count,
                                   &sendBuff[0],
                                   &recvCounts[0],
                                   &recvBuff[0]
                                 );

    // normalize
    double domainArea = grid.getGridDomainArea();
    double l2norm = sqrt ( recvBuff[0] ) / domainArea;

    return l2norm;
}
// =============================================================================
// Count the number of vortices by the total phase change along the boundary
// of the domain.
// TODO Make this work in multicore environments.
// Idea: Cauchy's integral formula: just caluculate the integral.
//       Numerically difficult when too close to origin (rather: points
//       closest to and furthest from origin too far apart).
int
Ginla::Helpers::
getVorticity ( const ComplexVector        & psi,
               const Recti::Grid::General & grid )
{
    int numProcs = psi.getMap()->getComm()->getSize();
    if ( numProcs!=1 )
        return -1;

    Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

    Teuchos::Array<int> boundaryIndices = grid.getBoundaryIndices();

    unsigned int n = boundaryIndices.length();

    // winding number algorithm
    // http://www.google.be/url?sa=t&source=web&ct=res&cd=1&ved=0CAkQFjAA&url=http%3A%2F%2Fwww.engr.colostate.edu%2F~dga%2Fdga%2Fpapers%2Fpoint_in_polygon.pdf&ei=ntxvS6eyOYLI-Qao2vjoAQ&usg=AFQjCNF_btLpRhTfkUQt34nXfsRylaF95g&sig2=Hh8I-SjpWC6queaj3QMeZw
    // TODO make vorticity INT
    double vorticity = 0;

    const double tol= 1.0e-15;

    double_complex currentZ;
    double_complex nextZ = psiView[ boundaryIndices[0] ];
    for ( unsigned int k=0; k<n; k++ )
    {
        currentZ = nextZ;
        if ( k<n-1 )
            nextZ = psiView[ boundaryIndices[k+1] ];
        else
            nextZ = psiView[ boundaryIndices[0] ];

        bool isK0Neg = std::imag ( currentZ ) < 0.0;
        bool isK1Neg = std::imag ( nextZ )    < 0.0;
        // check crossings of the negative real axis
        if ( isK0Neg ^ isK1Neg ) // crosses real axis
        {
            // calculate the intersection point
            double r = std::real ( currentZ ) + std::imag ( currentZ ) * ( std::real ( nextZ ) - std::real ( currentZ ) ) / ( std::imag ( currentZ )- std::imag ( nextZ ) );
            if ( std::signbit ( r ) ) // is negative?
            {
                if ( isK0Neg ) // crossing clockwise
                    --vorticity;
                else  // crossing counter-clockwise
                    ++vorticity;
            }
        }
        else if ( fabs ( std::imag ( currentZ ) ) < tol && std::real ( currentZ ) < 0.0 )
        {
            if ( isK1Neg )
                vorticity += 0.5;
            else
                vorticity -= 0.5;
        }
        else if ( fabs ( std::imag ( nextZ ) ) < tol && std::real ( nextZ ) < 0.0 )
        {
            if ( isK0Neg )
                vorticity -= 0.5;
            else
                vorticity += 0.5;
        }
    }

    return round ( vorticity );
}
// ============================================================================