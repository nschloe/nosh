/********************************************//**
 * The Ginzburg--Landau equations.
 ***********************************************/
#ifndef GINZBURGLANDAU_H
#define GINZBURGLANDAU_H

#include <complex>

#include "glBoundaryConditionsVirtual.h"
#include "GridVirtual.h"
#include "GridSquare.h"
#include "MagneticVectorPotential.h"

#include <Teuchos_ParameterList.hpp>

#include <Teuchos_XMLObject.hpp>

#include <Epetra_Map.h>

typedef std::complex<double> double_complex;
typedef Tpetra::Vector<double_complex,int> ComplexVector;

class GinzburgLandau
{
  public:

     /*! Default constructor.*/
     GinzburgLandau ( Teuchos::RCP<GridUniformVirtual>          &grid,
                      Teuchos::RCP<MagneticVectorPotential>     &potential,
                      Teuchos::RCP<GlBoundaryConditionsVirtual> &bc    );

     /*! Destructor. */
     ~GinzburgLandau();

     void
     setH0(const double h0);

     int
     getNumUnknowns() const;

     void
     setScaling( const double scaling);

     ComplexVector
     computeGlVector( const Teuchos::RCP<ComplexVector> & psi ) const;

     /*! Returns the coefficients of the jacobian system associated with the
         Ginzburg--Landau equations. */
     void
     getJacobianRow( const int                         eqnum,
                     const Teuchos::RCP<ComplexVector> &psi,
                     std::vector<int>                  &columnIndicesPsi,
                     std::vector<double_complex>       &valuesPsi,
                     std::vector<int>                  &columnIndicesPsiConj,
		     std::vector<double_complex>       &valuesPsiConj
                   ) const;

     /*! Get sparsity pattern of the jacobian system. */
     void
     getJacobianRowSparsity( int              eqnum,
                             std::vector<int> &columnIndicesPsi,
                             std::vector<int> &columnIndicesPsiConj
                           ) const;

     /*! Calcuate the grid approximation of the Gibbs free energy
       \f[
       \mathcal{G} = \int\nolimits_{\Omega} |\psi|^4 \,\mathrm{d}\omega
       \f]
       of a given state \f$\psi\f$. */
     double
     freeEnergy ( const ComplexVector &psi ) const;
     
     double
     normalizedScaledL2Norm ( const ComplexVector &psi ) const;

     /*! Calculate the vorticity of the current solution. */
     int
     getVorticity ( const ComplexVector &psi ) const;

     /*! Writes a solution \c psi to a file with all parameters that
      *  may be interesting. */
     void
     writeSolutionToFile( const Teuchos::RCP<const ComplexVector> &psi,
                          const std::string &filePath
                        ) const;

     /*! Writes an abstract state \c psi to a file (e.g., an eigenstate). The
      * only difference compared to the method \c writeSolutionToFile is
      * the number and the kind of parameters that are returned in the file.
      */
     void
     writeAbstractStateToFile( const Teuchos::RCP<const ComplexVector> &psi,
                               const std::string &filePath
                             ) const;

     /*! Appends useful statistics about a given state \c psi to the \c ofstream
      *  \c filestream.
      */
     void
     appendStats( std::ofstream & fileStream,
                  const bool header = false,
                  const Teuchos::RCP<const ComplexVector> &psi = Teuchos::null
                ) const;
      
  private:

     /*! Evaluates the Ginzburg--Landau equations.
         @param eqnum Index of the equation to evaluate.
         @param psi   Current order parameter \f$\psi\f$. */
     double_complex
     computeGl( const int           eqnum,
                const ComplexVector &psi
              ) const;
    
      //! Equation type enumerator.
      /*! Semantically separates the different types of conditions which must
          be applied at different parts of the rectangular grid. */
      enum equationType
      {
        BOUNDARY,
        INTERIOR,
        PHASE_CONDITION
      };

      const Teuchos::RCP<GridUniformVirtual>      grid_;
      const Teuchos::RCP<MagneticVectorPotential> A_;

      void getEquationType ( const int           eqnum,
                             equationType        &eqType,
                             int                 &eqIndex
                           ) const;
      
      /*! Have \c boundaryConditions_ declared as pointer as its class
          \c GlBoundaryConditionsVirtual is only available via a forward
	  declaration. */
      Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions_;

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

      void
      writeStateToFile( const Teuchos::RCP<const ComplexVector> &psi,
                        Teuchos::ParameterList &params,
                        const std::string &filePath) const;

};
#endif // GINZBURGLANDAU_H


void
readStateFromFile ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
		    const std::string           & filePath,
                    Teuchos::RCP<ComplexVector> & psi,
                    Teuchos::RCP<GridUniformVirtual>   & grid,
                    Teuchos::ParameterList      & params
                  );
