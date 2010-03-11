/********************************************//**
 * The Ginzburg--Landau equations.
 ***********************************************/
#ifndef GINZBURGLANDAU_H
#define GINZBURGLANDAU_H

#include "Ginla_Operator_Virtual.h"
#include "Recti_Grid_Abstract.h"

#include "Ginla_MagneticVectorPotential_Centered.h"

#include "Ginla_Perturbation_Virtual.h"

#include "Ginla_StatsWriter.h"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Epetra_Map.h>
#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

typedef std::complex<double> double_complex;
typedef Tpetra::Vector<double_complex,Thyra::Ordinal> ComplexVector;

class GinzburgLandau
{
public:

    /*! Default constructor.*/
    GinzburgLandau ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
                     const Teuchos::RCP<Ginla::StatsWriter>       & statsWriter,
                     const std::string                         & outputFormat = "VTI"
                   );
                   
    GinzburgLandau ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
                     const std::string                         & outputFormat = "VTI"
                   );

    /*! Constructor with a specified perturbation.*/
    GinzburgLandau ( const Teuchos::RCP<Ginla::Operator::Virtual>     & glOperator,
	             const Teuchos::RCP<Ginla::Perturbation::Virtual> & perturbation,
                     const std::string                                & outputFormat = "VTI"
                   );

    /*! Destructor. */
    ~GinzburgLandau();

    void 
    setParameters( const LOCA::ParameterVector & p );

    int
    getNumUnknowns() const;
    
    Teuchos::RCP<Ginla::StatsWriter>
    getStatsWriter();

    Teuchos::RCP<ComplexVector>
    computeGlVector ( const Teuchos::RCP<const ComplexVector> & psi ) const;

    /*! Returns the coefficients of the Jacobian system associated with the
        Ginzburg--Landau equations. */
    void
    getJacobianRow ( const int                           eqnum,
                     const Teuchos::RCP<ComplexVector> & psi,
                     Teuchos::Array<int>               & columnIndicesPsi,
                     Teuchos::Array<double_complex>    & valuesPsi,
                     Teuchos::Array<int>               & columnIndicesPsiConj,
                     Teuchos::Array<double_complex>    & valuesPsiConj
                   ) const;

    Teuchos::RCP<Ginla::Operator::Virtual>
    getOperator() const;

private:

    /*! Evaluates the Ginzburg--Landau equations.
        @param eqnum Index of the equation to evaluate.
        @param psi   Current order parameter \f$\psi\f$. */
    double_complex
    computeGl ( const int           eqnum,
                const ComplexVector &psi
              ) const;

    double_complex
    getInterior ( const int k,
                  const ComplexVector & psi
                ) const;

private:
    const Teuchos::RCP<Ginla::Operator::Virtual> glOperator_;

    const Teuchos::RCP<Ginla::Perturbation::Virtual> perturbation_;
    
    const Teuchos::RCP<Ginla::StatsWriter> statsWriter_;
    
    const std::string outputFormat_;

private:
    /*! Calculated the coefficients of the jacobian system associated with the
        Ginzburg--Landau equations. */
    void computeJacobianRow ( const bool                        fillValues,
                              const int                         eqnum,
                              const Teuchos::RCP<ComplexVector> &psi,
                              std::vector<int>                  &columnIndicesPsi,
                              std::vector<double_complex>       &valuesPsi,
                              std::vector<int>                  &columnIndicesPsiConj,
                              std::vector<double_complex>       &valuesPsiConj
                            ) const;

};
#endif // GINZBURGLANDAU_H
