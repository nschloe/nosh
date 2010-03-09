#include "ginzburgLandau.h"

#include "Ginla_Helpers.h"

// needed for reduceAllAndScatter in freeEnergy()
#include <Teuchos_Comm.hpp>



// =============================================================================
// Class constructor
GinzburgLandau::
GinzburgLandau ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
                 const Teuchos::RCP<Ginla::StatsWriter>       & statsWriter,
                 const std::string                         & outputFormat
               ) :
        glOperator_ ( glOperator ),
        perturbation_ ( Teuchos::null ),
        statsWriter_( statsWriter ),
        outputFormat_( outputFormat )
{
}
// =============================================================================
// Class constructor
GinzburgLandau::
GinzburgLandau ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
                 const std::string                         & outputFormat
               ) :
        glOperator_ ( glOperator ),
        perturbation_ ( Teuchos::null ),
        statsWriter_( Teuchos::null ),
        outputFormat_( outputFormat )
{
}
// =============================================================================
// Class constructor containing a perturbation
GinzburgLandau::
GinzburgLandau ( const Teuchos::RCP<Ginla::Operator::Virtual>     & glOperator,
                 const Teuchos::RCP<Ginla::Perturbation::Virtual> & perturbation,
                 const std::string                             & outputFormat
               ) :
        glOperator_ ( glOperator ),
        perturbation_ ( perturbation ),
        outputFormat_( outputFormat )
{
}
// =============================================================================
// Destructor
GinzburgLandau::~GinzburgLandau()
{
}
// =============================================================================
void
GinzburgLandau::setParameters ( const LOCA::ParameterVector & p )
{
    TEUCHOS_ASSERT( glOperator_.is_valid_ptr() && !glOperator_.is_null() );
    glOperator_->setParameters ( p );
    
    if ( perturbation_.is_valid_ptr() && !perturbation_.is_null() )
        perturbation_->setParameters ( p );
}
// =============================================================================
int
GinzburgLandau::getNumUnknowns() const
{
    return glOperator_->getGrid()->getNumGridPoints();
}
// =============================================================================
Teuchos::RCP<Ginla::StatsWriter>
GinzburgLandau::getStatsWriter()
{
    return statsWriter_;
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

        if ( !perturbation_.is_null() )
            glVecView[k] += perturbation_->computePerturbation ( globalIndex );
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
    glOperator_->updatePsi ( psi );
    glOperator_->getJacobianRow ( eqnum,
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
// Idea: Cauchy's integral formula: just caluculate the integral.
//       Numerically difficult when too close to origin (rather: points
//       closest to and furthest from origin too far apart).
int
GinzburgLandau::getVorticity ( const ComplexVector &psi
                             ) const
{
    int numProcs = psi.getMap()->getComm()->getSize();
    if ( numProcs!=1 )
        return -1;

    Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

    Teuchos::Array<int> boundaryIndices =
        glOperator_->getGrid()->getBoundaryIndices();

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
// =============================================================================
void
GinzburgLandau::
writeStateToFile ( const Teuchos::RCP<const ComplexVector> & psi,
                   LOCA::ParameterVector                   & params,
                   const std::string                       & fileBaseName
                 ) const
{
    Teuchos::RCP<Teuchos::ParameterList> p =
        Ginla::Helpers::locaParameterVector2teuchosParameterList( params );
        
    std::string filenameExtension;
    if ( outputFormat_.compare("VTI")==0 )
      filenameExtension = "vti";
    else if ( outputFormat_.compare("VTK")==0 )
      filenameExtension = "vtk";
    else
      TEST_FOR_EXCEPTION( true,
                          std::runtime_error,
                          "outputFormat_ (\"" << outputFormat_
                          << "\") must be either one of \"VTI\", \"VTK\"" );
        
    std::string fileName = fileBaseName + "." + filenameExtension;
    glOperator_->getGrid()->writeWithGrid ( *psi, *p, fileName );
    return;
}
// =============================================================================
void
GinzburgLandau::
writeSolutionToFile ( const Teuchos::RCP<const ComplexVector> & psi,
                      const std::string                       & fileBaseName
                    ) const
{
    // create a parameter list that contains useful items for a solution file
    Teuchos::RCP<LOCA::ParameterVector> params = glOperator_->getParameters();
    writeStateToFile ( psi, *params, fileBaseName );
    return;
}
// =============================================================================
void
GinzburgLandau::
writeAbstractStateToFile ( const Teuchos::RCP<const ComplexVector> & psi,
                           const std::string                       & fileBaseName
                         ) const
{
    LOCA::ParameterVector params;
    writeStateToFile ( psi, params, fileBaseName );
    return;
}
// =============================================================================
void
GinzburgLandau::
appendStats ( const Teuchos::RCP<const ComplexVector> & psi
            ) const
{ 
    TEUCHOS_ASSERT( statsWriter_.is_valid_ptr() && !statsWriter_.is_null() );
  
    // put the parameter list into statsWriter_
    std::string labelPrepend = "1";
    Ginla::Helpers::appendToTeuchosParameterList( *(statsWriter_->getList()),
                                               *(glOperator_->getParameters()),
                                               labelPrepend );
    
    TEUCHOS_ASSERT( psi.is_valid_ptr() && !psi.is_null() );
    
    statsWriter_->getList()->set( "2free energy", freeEnergy ( *psi ) );
    statsWriter_->getList()->set( "2||x||_2 scaled", normalizedScaledL2Norm ( *psi ) );
    statsWriter_->getList()->set( "2vorticity", getVorticity ( *psi ) );

    return;
}
// =============================================================================
