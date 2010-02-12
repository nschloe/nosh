#include "ginzburgLandau.h"

#include <Teuchos_RCP.hpp>

// really needed?
// --> reduceAllAndScatter in freeEnergy()
#include <Teuchos_Comm.hpp>

#include <Tpetra_Map.hpp>

// complex unit
const double_complex I ( 0,1 );

// =============================================================================
// Class constructor
GinzburgLandau::GinzburgLandau ( const Teuchos::RCP<GlOperatorVirtual>  & glOperator
                               ) :
        glOperator_ ( glOperator )
{
}
// =============================================================================
// Destructor
GinzburgLandau::~GinzburgLandau()
{
}
// =============================================================================
void
GinzburgLandau::setH0 ( const double h0 )
{
    glOperator_->setH0 ( h0 );
}
// =============================================================================
// TODO delete?
double
GinzburgLandau::getH0() const
{
    return glOperator_->getH0();
}
// =============================================================================
void
GinzburgLandau::setScaling ( const double scaling )
{
    glOperator_->setScaling ( scaling );
}
// =============================================================================
void
GinzburgLandau::setChi ( const double chi )
{
    glOperator_->setChi ( chi );
}
// =============================================================================
int
GinzburgLandau::getNumUnknowns() const
{
    return glOperator_->getGrid()->getNumGridPoints();
}
// =============================================================================
Teuchos::RCP<ComplexVector>
GinzburgLandau::computeGlVector ( const Teuchos::RCP<const ComplexVector> & psi
                                ) const
{
    TEST_FOR_EXCEPTION ( !psi.is_valid_ptr() || psi.is_null(),
                         std::logic_error,
                         "Input argument 'psi' not properly initialized." );

    // setup output vector with the same map as psi
    Teuchos::RCP<ComplexVector> glVec =
        Teuchos::rcp ( new ComplexVector ( psi->getMap(), true ) );

    glOperator_->updatePsi ( psi );

    Teuchos::ArrayRCP<double_complex> glVecView = glVec->get1dViewNonConst();
    // loop over the nodes
    for ( unsigned int k=0; k<psi->getLocalLength(); k++ )
    {
        int globalIndex = psi->getMap()->getGlobalElement ( k );
        glVecView[k] = glOperator_->getEntry ( globalIndex );
    }

    return glVec;
}
// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
void
GinzburgLandau::getJacobianRow ( const int                           eqnum,
                                 const Teuchos::RCP<ComplexVector> & psi,
                                 Teuchos::Array<int>               & columnIndicesPsi,
                                 Teuchos::Array<double_complex>    & valuesPsi,
                                 Teuchos::Array<int>               & columnIndicesPsiConj,
                                 Teuchos::Array<double_complex>    & valuesPsiConj
                               ) const
{
    glOperator_->updatePsi( psi );
    glOperator_->getJacobianRow( eqnum,
                                 columnIndicesPsi, valuesPsi,
                                 columnIndicesPsiConj, valuesPsiConj );
    return;
}
// =============================================================================
// calculate the free energy of a state
double
GinzburgLandau::freeEnergy ( const ComplexVector & psi
                           ) const
{
    double localEnergy = 0.0;

    Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

    // sum up the energy on each processor
    for ( unsigned int k=0; k<psi.getLocalLength(); k++ )
    {
        int kGlobal = psi.getMap()->getGlobalElement ( k );
        double area = glOperator_->getGrid()->cellArea ( kGlobal );
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
    double maxEnergy = 1.0 * glOperator_->getGrid()->getGridDomainArea();
    double globalEnergy = recvBuff[0] / maxEnergy;

    return globalEnergy;
}
// =============================================================================
double
GinzburgLandau::normalizedScaledL2Norm ( const ComplexVector &psi
                                       ) const
{
    double localSum = 0.0;

    Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

    // sum up on each processor
    for ( unsigned int k=0; k<psi.getLocalLength(); k++ )
    {
        int kGlobal = psi.getMap()->getGlobalElement ( k );
        double area = glOperator_->getGrid()->cellArea ( kGlobal );
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
    double domainArea = glOperator_->getGrid()->getGridDomainArea();
    double l2norm = sqrt ( recvBuff[0] ) / domainArea;

    return l2norm;
}
// =============================================================================
// Count the number of vortices by the total phase change along the boundary
// of the domain.
// TODO Make this work in multicore environments.
int
GinzburgLandau::getVorticity ( const ComplexVector &psi
                             ) const
{
    return 0;
//     int numProcs = psi.getMap()->getComm()->getSize();
//     if ( numProcs!=1 )
//         return -1;
//
//     int vorticity = 0;
//
//     const double PI = 3.14159265358979323846264338327950288419716939937510;
//     // Consider jumps in the argument greater than threshold phase jumps.
//     const double threshold = 1.5*PI;
//
//     // Get a view of the whole vector.
//     // Remember: This only works with one core.
//     Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();
//
//     int numBoundaryPoints = grid_->getNumBoundaryPoints();
//
//     int psiLowerOffset = psiView.lowerOffset();
//     int psiUpperOffset = psiView.upperOffset();
//     int boundaryIndex = 0;
//     int globalIndex = grid_->boundaryIndex2globalIndex ( boundaryIndex );
//
//     TEST_FOR_EXCEPTION ( globalIndex < psiLowerOffset || globalIndex > psiUpperOffset,
//                          std::out_of_range,
//                          "Index globalIndex=" << globalIndex << " out of bounds (lower="
//                          << psiLowerOffset << ", upper=" << psiUpperOffset << ")" );
//
//     double angle = std::arg ( psiView[globalIndex] );
//     double angle0 = angle;
//     double anglePrevious;
//     for ( boundaryIndex=1; boundaryIndex<numBoundaryPoints; boundaryIndex++ )
//     {
//         anglePrevious = angle;
//         globalIndex   = grid_->boundaryIndex2globalIndex ( boundaryIndex );
//         TEST_FOR_EXCEPTION ( globalIndex < psiLowerOffset || globalIndex > psiUpperOffset,
//                              std::out_of_range,
//                              "Index globalIndex=" << globalIndex << " out of bounds (lower="
//                              << psiLowerOffset << ", upper=" << psiUpperOffset << ")" );
//         angle = std::arg ( psiView[globalIndex] );
//         if ( angle-anglePrevious<-threshold )
//             vorticity++;
//         else if ( angle-anglePrevious>threshold )
//             vorticity--;
//     }
//
//     // close the circle
//     anglePrevious = angle;
//     if ( angle0-anglePrevious<-threshold )
//         vorticity++;
//     else if ( angle0-anglePrevious>threshold )
//         vorticity--;
//
//     return vorticity;
}
// =============================================================================
void
GinzburgLandau::writeStateToFile ( const Teuchos::RCP<const ComplexVector> &psi,
                                   Teuchos::ParameterList &params,
                                   const std::string &filePath
                                 ) const
{
    glOperator_->getGrid()->writeWithGrid ( *psi, params, filePath );
}
// =============================================================================
void
GinzburgLandau::writeSolutionToFile ( const Teuchos::RCP<const ComplexVector> &psi,
                                      const std::string &filePath
                                    ) const
{
    // create a parameter list that contains useful items for a solution file
    Teuchos::ParameterList params;

    // TODO avoid calculating free energy and vorticity twice
    double energy = freeEnergy ( *psi );
    int vorticity = getVorticity ( *psi );

    params.get ( "H0", glOperator_->getH0() );
    params.get ( "free energy", energy );
    params.get ( "vorticity", vorticity );

    writeStateToFile ( psi, params, filePath );
}
// =============================================================================
void
GinzburgLandau::writeAbstractStateToFile ( const Teuchos::RCP<const ComplexVector> &psi,
        const std::string &filePath
                                         ) const
{
    Teuchos::ParameterList params;
    writeStateToFile ( psi, params, filePath );
}
// =============================================================================
void
GinzburgLandau::appendStats ( std::ofstream & fileStream,
                              const bool header,
                              const Teuchos::RCP<const ComplexVector> & psi
                            ) const
{
    if ( header )
    {
        fileStream
        << "H0                \t"
        << "scaling           \t"
        << "free energy       \t"
        << "||x||_2 scaled    \t"
        << "vorticity";
    }
    else
    {
        // TODO avoid calculating the free energy twice
        double energy = freeEnergy ( *psi );
        double l2norm = normalizedScaledL2Norm ( *psi );
        int vorticity = getVorticity ( *psi );
        fileStream
        << glOperator_->getH0() << " \t"
        << glOperator_->getScaling() << " \t"
        << energy << " \t"
        << l2norm << " \t"
        << vorticity;
    }

}
// =============================================================================
