/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl"omer

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

#include "Ginla_State.h"

// =============================================================================
Ginla::State::
State( const Teuchos::RCP<ComplexVector>              & psi,
       const Teuchos::RCP<const Recti::Grid::General> & grid ):
       psi_( psi ),
       chi_( 0.0 ),
       grid_( grid )
{
    TEUCHOS_ASSERT_EQUALITY( psi->getGlobalLength(),
                             grid->getNumGridPoints() );
}
// =============================================================================
Ginla::State::
State( const Teuchos::RCP<const ComplexMap>           & map,
       const Teuchos::RCP<const Recti::Grid::General> & grid ):
       chi_( 0.0 ),
       grid_( grid )
{
   bool zeroOut = true;
   psi_ = Teuchos::rcp( new ComplexVector( map, zeroOut ) );
}
// =============================================================================
Teuchos::RCP<const ComplexVector>
Ginla::State::
getValuesConst () const
{
    return psi_; 
}
// =============================================================================
Teuchos::RCP<ComplexVector>
Ginla::State::
getValuesNonConst ()
{
    return psi_; 
}
// =============================================================================
void
Ginla::State::
setValues ( Teuchos::RCP<ComplexVector> & psi )
{
    TEUCHOS_ASSERT_EQUALITY( psi->getGlobalLength(),
                             grid_->getNumGridPoints() );
    psi_ = psi; 
}
// =============================================================================
const Teuchos::RCP<const Recti::Grid::General>
Ginla::State::
getGrid () const
{
    return grid_; 
}
// =============================================================================
double
Ginla::State::
getChi () const
{
  return chi_;
}
// =============================================================================
double
Ginla::State::
freeEnergy () const
{
    double localEnergy = 0.0;

    Teuchos::ArrayRCP<const double_complex> psiView = psi_->get1dView();

    // sum up the energy on each processor
    for ( unsigned int k=0; k<psi_->getLocalLength(); k++ )
    {
        int kGlobal = psi_->getMap()->getGlobalElement ( k );
        double area = grid_->cellArea ( kGlobal );
        localEnergy -= area * pow ( norm ( psiView[k] ),2 );
    }

    // reduce and scatter such that energy is available on
    // all cores
    int count = 1; // send *one* integer
    int numProcs =  psi_->getMap()->getComm()->getSize();
    Teuchos::Array<double> sendBuff ( count ), recvBuff ( count );

    // fill send buffer
    sendBuff[0] = localEnergy;

    Teuchos::Array<int> recvCounts ( numProcs );
// fill recvCounts with {1,...,1}
    int numItemsPerProcess = 1;
    std::fill ( recvCounts.begin(), recvCounts.end(), numItemsPerProcess );

    Teuchos::reduceAllAndScatter ( * ( psi_->getMap()->getComm() ),
                                   Teuchos::REDUCE_SUM,
                                   count,
                                   &sendBuff[0],
                                   &recvCounts[0],
                                   &recvBuff[0]
                                 );

    // normalize
    double maxEnergy = 1.0 * grid_->getGridDomainArea();
    double globalEnergy = recvBuff[0] / maxEnergy;

    return globalEnergy;
}
// =============================================================================
double
Ginla::State::
normalizedScaledL2Norm () const
{
    double localSum = 0.0;

    Teuchos::ArrayRCP<const double_complex> psiView = psi_->get1dView();

    // sum up on each processor
    for ( unsigned int k=0; k<psi_->getLocalLength(); k++ )
    {
        int kGlobal = psi_->getMap()->getGlobalElement ( k );
        double area = grid_->cellArea ( kGlobal );
        localSum += area * norm ( psiView[k] );
    }

    // reduce and scatter such that energy is available on
    // all cores
    int count = 1; // send *one* integer
    Teuchos::Array<double> sendBuff ( count ), recvBuff ( count );

    // fill send buffer
    sendBuff[0] = localSum;

    int numProcs =  psi_->getMap()->getComm()->getSize();
    Teuchos::Array<int> recvCounts ( numProcs );
// fill recvCounts with {1,...,1}
    int numItemsPerProcess = 1;
    std::fill ( recvCounts.begin(), recvCounts.end(), numItemsPerProcess );

    Teuchos::reduceAllAndScatter ( * ( psi_->getMap()->getComm() ),
                                   Teuchos::REDUCE_SUM,
                                   count,
                                   &sendBuff[0],
                                   &recvCounts[0],
                                   &recvBuff[0]
                                 );

    // normalize
    double domainArea = grid_->getGridDomainArea();
    double l2norm = sqrt ( recvBuff[0] ) / domainArea;

    return l2norm;
}
// =============================================================================
void
Ginla::State::
update( const double alpha,
        const Ginla::State b,
        const double beta  )
{
  psi_->update( alpha, *(b.getValuesConst()), beta );
  chi_ = alpha*chi_ + beta * b.getChi();
}
// =============================================================================
// Count the number of vortices by the total phase change along the boundary
// of the domain.
// TODO Make this work in multicore environments.
// Idea: Cauchy's integral formula: just caluculate the integral.
//       Numerically difficult when too close to origin (rather: points
//       closest to and furthest from origin too far apart).
int
Ginla::State::
getVorticity () const
{
    int numProcs = psi_->getMap()->getComm()->getSize();
    if ( numProcs!=1 )
        return -1;

    Teuchos::ArrayRCP<const double_complex> psiView = psi_->get1dView();

    Teuchos::Array<int> boundaryIndices = grid_->getBoundaryIndices();

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