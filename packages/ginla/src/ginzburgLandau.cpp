#include "ginzburgLandau.h"

#include "Ginla_Helpers.h"

// needed for reduceAllAndScatter in freeEnergy()
#include <Teuchos_Comm.hpp>


// =============================================================================
// Class constructor
GinzburgLandau::
GinzburgLandau ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
                 const Teuchos::RCP<Ginla::StatsWriter>       & statsWriter,
                 const std::string                            & outputFormat
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
                 const std::string                            & outputFormat
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
    
    statsWriter_->getList()->set( "2free energy", Ginla::Helpers::freeEnergy ( *psi, *(glOperator_->getGrid()) ) );
    statsWriter_->getList()->set( "2||x||_2 scaled", Ginla::Helpers::normalizedScaledL2Norm ( *psi, *(glOperator_->getGrid()) ) );
    statsWriter_->getList()->set( "2vorticity", Ginla::Helpers::getVorticity ( *psi, *(glOperator_->getGrid()) ) );

    return;
}
// =============================================================================
