/*
 * GlSystemWithConstraint.h
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl\"omer
 */
#ifndef GINLA_LOCASYSTEM_BORDERED_H
#define GINLA_LOCASYSTEM_BORDERED_H

#include "Ginla_LocaSystem_Virtual.h"

#include <Epetra_Import.h>

#include "Ginla_LocaSystem_Default.h"



// forward declarations
namespace Ginla {
  namespace Operator {
    class Virtual;
  }
  namespace IO {
    class StatsWriter;
    class StateWriter;
  }
}
namespace LOCA {
  class Stepper;
}


namespace Ginla {
  namespace LocaSystem {

class Bordered:
            public Ginla::LocaSystem::Virtual
{
public:

    //! Constructor with initial guess.
    Bordered ( const Teuchos::RCP<Ginla::Operator::Virtual>           & glOperator,
               const Teuchos::RCP<const Epetra_Comm>                  & eComm,
               const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & complexMap,
               const Teuchos::RCP<Ginla::IO::StatsWriter>             & statsWriter,
               const Teuchos::RCP<Ginla::IO::StateWriter>             & stateWriter );

    //! Destructor
    ~Bordered();

    //! Evaluate the Ginzburg--Landau functions at a given state defined
    //! by the input vector x.
    virtual bool
    computeF ( const Epetra_Vector &x,
               Epetra_Vector &F,
               const NOX::Epetra::Interface::Required::FillType fillFlag = Residual );

    //! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
    //! at a given state defined by the input vector x.
    virtual bool
    computeJacobian ( const Epetra_Vector &x, Epetra_Operator &Jac );

    virtual bool
    computeShiftedMatrix ( double alpha,
                           double beta,
                           const Epetra_Vector &x,
                           Epetra_Operator &A );

    //! Dummy preconditioner function. So far does nothing but throwing
    //! an exception when called.
    virtual bool
    computePreconditioner ( const Epetra_Vector    & x,
                            Epetra_Operator        & Prec,
                            Teuchos::ParameterList * precParams = 0 );
    
    //! Returns the current Jacobian.
    //! @return Reference-counted pointer to the Jacobian.
    Teuchos::RCP<Epetra_CrsMatrix>
    getJacobian() const;

    Teuchos::RCP<Epetra_CrsMatrix>
    getPreconditioner() const;

    //! Set the problem parameters.
    virtual void
    setParameters ( const LOCA::ParameterVector & p );

    //! Set directory to where all output gets written.
    //! @param directory Name of the directory.
    void
    setOutputDir ( const string & directory );

    //! Print the solution x along with the continuation parameter conParam
    //! to a file. This function is called internally by LOCA to print
    //! solutions of each continuation step.
    virtual void
    printSolution ( const Epetra_Vector & x,
                    double conParam );
        
    //! Used to print eigenvectors.
    virtual void
    printSolution ( const Epetra_Vector & x,
                    const std::string   & filenameAppendix
                  ) const;

    virtual void
    setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper );

    // This function is necessary to break the circular dependency with the
    // LOCA_Stepper object to allow for a clean termination
    virtual void
    releaseLocaStepper();

    virtual Teuchos::RCP<const Komplex>
    getKomplex() const;

    virtual Teuchos::RCP<const Epetra_Map>
    getMap() const;

    //! Creates a state suitable for usage in the present interface.
    //! This function can be used, for example, to create an initial guess for
    //! the system.
    virtual Teuchos::RCP<Epetra_Vector>
    createSystemVector( const Teuchos::ParameterList & p );
                          
    //! Extracts a complex-valued state from a vector of the linear
    //! equation system.
    //! Can be used to tranform eigenvalues into states from within
    //! other classes.
//     Teuchos::RCP<ComplexVector>
//     extractPsi( const Epetra_Vector & x ) const;

    //! Creates a valid state from a raw real-valued EpetraVector.
    Teuchos::RCP<Ginla::State>
    createState( const Epetra_Vector & x ) const;

private:

    void
    fillBorderedMatrix ( const Teuchos::RCP<      Epetra_CrsMatrix> & extendedMatrix,
                         const Teuchos::RCP<const Epetra_CrsMatrix> & regularMatrix,
                         const Epetra_Vector                        & rightBorder,
                         Epetra_Vector                              & lowerBorder,
                         double                                       d,
                         bool                                         firstTime
                       ) const;

    int
    PutRow ( const Teuchos::RCP<Epetra_CrsMatrix> A,
             const int      Row,
             const int      numIndices,
             double       * values,
             int          * indices,
             const bool     firstTime ) const;


    //! Creates a map identical to the input argument \c realMap, but
    //! extended by one entry.
    //! This extra slot is usually used for the phase condition.
    Teuchos::RCP<Epetra_Map>
    createExtendedRealMap ( const Epetra_BlockMap & realMap ) const;

    void
    createJacobian ( const Epetra_Vector &x );

    //! Print method for the continuation in one parameter.
    void
    printSolutionOneParameterContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                          ) const;

    //! Print method for turning point continuation continuation.
    void
    printSolutionTurningPointContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                          ) const;


    //! Write statistics about the current continuation step to the file
    //! \c outputDataFileName_ .
    void
    writeContinuationStats ( const int conStep,
                             const Teuchos::RCP<const ComplexVector> psi ) const;

    //! Translate an Epetra_Comm into a Teuchos::Comm<int>, no matter the Thyra::Ordinal.
    Teuchos::RCP<const Teuchos::Comm<int> >
    create_CommInt ( const Teuchos::RCP<const Epetra_Comm> &epetraComm );

private:

    Ginla::LocaSystem::Default glSystem_;

    Teuchos::RCP<const Epetra_BlockMap> regularMap_;
    Teuchos::RCP<Epetra_Map> extendedMap_;
    const Teuchos::RCP<Epetra_CrsMatrix> jacobian_;
    Teuchos::RCP<Epetra_CrsMatrix> preconditioner_;

    std::string stepNumFileNameFormat_;

    bool firstTime_;

    const Epetra_Import importFromExtendedMap_;
    const Epetra_Import importFromRegularMap_;
};

  } // namespace LocaSystem
} // namespace GL

#endif // GINLA_LOCASYSTEM_BORDERED_H
