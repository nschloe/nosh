/********************************************//**
 * The Ginzburg--Landau equations.
 ***********************************************/
#include <complex>
#include <string> // for the file output functions

#include "staggeredGrid.h"

#include <Epetra_Map.h>

#include <Teuchos_ParameterList.hpp>

#include <Teuchos_XMLObject.hpp>

// abbreviate the complex type name
typedef std::complex<double> double_complex;

class GinzburgLandau
{
  public:

     /*! Default constructor.
         @param nx Number of boxes in each spatial dimension.
         @param edgelength Edge length of the square domain.
         @param h0 (Initial) external magnetic field strength. */
     GinzburgLandau( int nx,
                     double edgelength,
                     double h0 );

     /*! Destructor. */
     ~GinzburgLandau();

     /*! Returns a pointer to the \f$A\f$ grid in use.*/
     StaggeredGrid::StaggeredGrid* getStaggeredGrid();

     /*! Evaluates the Ginzburg--Landau equations.
         @param eqnum Index of the equation to evaluate.
         @param psi   Current order parameter \f$\psi\f$. */
     double_complex computeGl( const int                         eqnum,
                               const std::vector<double_complex> &psi   );

     /*! Returns the coefficients of the jacobian system associated with the
         Ginzburg--Landau equations. */
     void getJacobianRow( const int                         eqnum,
                          const std::vector<double_complex> &psi,
                          std::vector<int>                  &columnIndicesPsi,
                          std::vector<double_complex>       &valuesPsi,
                          std::vector<int>                  &columnIndicesPsiConj,
                          std::vector<double_complex>       &valuesPsiConj );

     /*! Get sparsity pattern of the jacobian system. */
     void getJacobianRowSparsity( int              eqnum,
                                  std::vector<int> &columnIndicesPsi,
                                  std::vector<int> &columnIndicesPsiConj );

     /*! Calcuate the grid approximation of the Gibbs free energy
       \f[
       \mathcal{G} = \int\nolimits_{\Omega} |\psi|^4 \,\mathrm{d}\omega
       \f]
       of a given state \f$\psi\f$. */
     double freeEnergy( const std::vector<double_complex> &psi );

     /*! Count the number of vortices. */
     int countVortices ( const std::vector<double_complex> &psi );

     /*! Print the solution \f$\psi\f$ to a legacy VTK file for viewing
         with ParaView, for example. */
     void psiToLegacyVtkFile( const std::vector<double_complex> &psi,
                              const std::string                 &filename );

     /*! Print the solution \f$\psi\f$ to an (XML-style) VTK file for viewing
         with ParaView, for example. */
     void psiToVtkFile( const std::vector<double_complex> &psi,
                        const Teuchos::ParameterList      &problemParams,
                        const std::string                 &filename       );

     void vtkFileToPsi( const std::string           &filename );
//                         std::vector<double_complex> *psi,
//                         Teuchos::ParameterList      *problemParams

     void psiToXdmfFile( const std::vector<double_complex> &psi,
                         const std::string                 &filename,
                         const Epetra_Map                  &StandardMap,
                         const Epetra_Comm                 &comm         );

  private:

      const Teuchos::XMLObject* xmlFind ( const Teuchos::XMLObject *xmlObj,
                                            const std::string        tag      );

      const Teuchos::XMLObject* xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                                                   const std::string        tag,
                                                   const std::string        attribute,
                                                   const std::string        value      );

      StaggeredGrid::StaggeredGrid sGrid;

      //! Equation type enumerator.
      /*! Semantically separates the different types of conditions which must
          be applied at different parts of the rectangular grid. */
      enum equationType
      {
        BOTTOMLEFT,
        BOTTOMRIGHT,
        TOPLEFT,
        TOPRIGHT,
        BOTTOM,
        TOP,
        LEFT,
        RIGHT,
        INTERIOR,
        PHASE_CONDITION
      };

      void getEquationType( const int,
                            equationType&,
                            int* );

      int getKLeft ( const int* i );
      int getKRight( const int* i );
      int getKBelow( const int* i );
      int getKAbove( const int* i );

      enum filltype { VALUES, SPARSITY };

      /*! Calculated the coefficients of the jacobian system associated with the
          Ginzburg--Landau equations. */
      void computeJacobianRow( const filltype                    ft,
                               const int                         eqnum,
                               const std::vector<double_complex> &psi,
                               std::vector<int>                  &columnIndicesPsi,
                               std::vector<double_complex>       &valuesPsi,
                               std::vector<int>                  &columnIndicesPsiConj,
                               std::vector<double_complex>       &valuesPsiConj         );

};