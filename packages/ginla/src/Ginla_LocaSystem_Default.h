/*! Jacobian system for the Ginzburg--Landau problem.
 *  This routine can be used as an interface to NOX.
 ******************************************************************************/
#ifndef GLSYSTEM_H
#define GLSYSTEM_H

#include "Ginla_Komplex.h"
#include "Ginla_IO_StatsWriter.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_Operator_Virtual.h"
#include "Ginla_Perturbation_Virtual.h"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>

#include <Epetra_CrsMatrix.h>

#include <Teuchos_ParameterList.hpp>

#include <NOX_Epetra_Interface_Required.H> // NOX base class
#include <NOX_Epetra_Interface_Jacobian.H> // NOX base class
#include <NOX_Epetra_Interface_Preconditioner.H> // NOX base class
#include <LOCA_Epetra_Interface_Required.H> // LOCA base class
#include <LOCA_Epetra_Interface_TimeDependent.H> // LOCA base class
#include <NOX_Abstract_PrePostOperator.H>

#include <LOCA_Parameter_Vector.H>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

#include <NOX_Abstract_Group.H>

#include <LOCA_Stepper.H>

typedef Tpetra::Vector<double_complex, Thyra::Ordinal> ComplexVector;

namespace Ginla {
  namespace LocaSystem {

class Default:
            public NOX::Epetra::Interface::Jacobian,
            public NOX::Epetra::Interface::Preconditioner,
            public LOCA::Epetra::Interface::TimeDependent
{
public:

    //! Constructor with initial guess.
    Default ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
              const Teuchos::RCP<const Epetra_Comm>        & eComm,
              const Teuchos::RCP<const ComplexVector>      & psi,
              const Teuchos::RCP<Ginla::IO::StatsWriter>   & statsWriter,
              const Teuchos::RCP<Ginla::IO::StateWriter>   & stateWriter
            );

    //! Destructor
    ~Default();

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
    computePreconditioner ( const Epetra_Vector &x, Epetra_Operator &Prec,
                            Teuchos::ParameterList* precParams = 0 );

    //! Returns the current state. Not necessarily a solution to the problem.
    //! @return Reference-counted pointer to the current state.
    Teuchos::RCP<Epetra_Vector>
    getSolution() const;

    //! Returns the current Jacobian.
    //! @return Reference-counted pointer to the Jacobian.
    Teuchos::RCP<Epetra_CrsMatrix>
    getJacobian() const;

    Teuchos::RCP<Epetra_CrsMatrix>
    getPreconditioner() const;

    //! Set the problem parameters.
    virtual void
    setParameters ( const LOCA::ParameterVector &p );

    //! Set directory to where all output gets written.
    //! @param directory Name of the directory.
    void
    setOutputDir ( const string &directory );

    //! Print the solution x along with the continuation parameter conParam
    //! to a file. This function is called internally by LOCA to print
    //! solutions of each continuation step.
    virtual void
    printSolution ( const Epetra_Vector &x, double conParam );

    void
    setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper );

    // This function is necessary to break the circular dependency with the
    // LOCA_Stepper object to allow for a clean termination
    void
    releaseLocaStepper();

    Teuchos::RCP<const Epetra_Map>
    getRealMap() const;

    Teuchos::RCP<const Ginla::Komplex>
    getGlKomplex() const;

private:

    enum continuationType
    {
        ONEPARAMETER,
        TURNINGPOINT
    };

    int
    numDigits ( const int i );

    //! Print method for the continuation in one parameter.
    void
    printSolutionOneParameterContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                          );

    //! Print method for turning point continuation continuation.
    void
    printSolutionTurningPointContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                          );

    //! Write statistics about the current continuation step to the file
    //! \c outputDataFileName_ .
    void
    writeContinuationStats ( const Teuchos::RCP<const ComplexVector> & psi );

    //! Translate an Epetra_Comm into a Teuchos::Comm<int>,
    //! no matter the Thyra::Ordinal.
    Teuchos::RCP<const Teuchos::Comm<int> >
    create_CommInt ( const Teuchos::RCP<const Epetra_Comm> &epetraComm );


private:
    const Teuchos::RCP<Ginla::Operator::Virtual> glOperator_;
    const Teuchos::RCP<Ginla::Perturbation::Virtual> perturbation_;
    Teuchos::RCP<Epetra_CrsMatrix> preconditioner_;

    Teuchos::RCP<const LOCA::Stepper> stepper_;

    Teuchos::RCP<Ginla::Komplex> glKomplex_;
    
    Teuchos::RCP<Epetra_Vector> initialSolution_;
    
    const Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter_;
    const Teuchos::RCP<Ginla::IO::StateWriter> stateWriter_;
    
    bool firstTime_;
    
    unsigned int maxNumDigits_; 
    
    continuationType continuationType_;
};

  } // namespace LocaSystem
} // namespace GL

#endif // GLSYSTEM_H
