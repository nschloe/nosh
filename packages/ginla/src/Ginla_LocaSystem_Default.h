/*! Jacobian system for the Ginzburg--Landau problem.
 *  This routine can be used as an interface to NOX.
 ******************************************************************************/
#ifndef GINLA_LOCASYSTEM_DEFAULT_H
#define GINLA_LOCASYSTEM_DEFAULT_H

#include "Ginla_LocaSystem_Virtual.h"

// forward declarations
namespace Ginla {
  namespace Operator {
    class Virtual;
  }
  namespace IO {
    class StatsWriter;
    class StateWriter;
  }
  namespace Perturbation {
    class Virtual;
  }
}
namespace LOCA {
  class Stepper;
}


namespace Ginla {
  namespace LocaSystem {

class Default:
            public Ginla::LocaSystem::Virtual
{
public:

    //! Constructor with initial guess.
    Default ( const Teuchos::RCP<Ginla::Operator::Virtual>           & glOperator,
              const Teuchos::RCP<const Epetra_Comm>                  & eComm,
              const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & complexMap,
              const Teuchos::RCP<Ginla::IO::StatsWriter>             & statsWriter,
              const Teuchos::RCP<Ginla::IO::StateWriter>             & stateWriter
            );

    //! Destructor
    ~Default();

    //! Evaluate the Ginzburg--Landau functions at a given state defined
    //! by the input vector x.
    virtual bool
    computeF ( const Epetra_Vector & x,
               Epetra_Vector       & F,
               const NOX::Epetra::Interface::Required::FillType fillFlag = Residual );
               
    //! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
    //! at a given state defined by the input vector x.
    virtual bool
    computeJacobian ( const Epetra_Vector & x,
                      Epetra_Operator     & Jac );

    virtual bool
    computeShiftedMatrix ( double alpha,
                           double beta,
                           const Epetra_Vector & x,
                           Epetra_Operator     & A );

    //! Dummy preconditioner function. So far does nothing but throwing
    //! an exception when called.
    virtual bool
    computePreconditioner ( const Epetra_Vector & x,
                            Epetra_Operator     & Prec,
                            Teuchos::ParameterList* precParams = 0 );

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
    printSolution ( const Epetra_Vector & x,
                    double conParam );
                    
    //! Used to print eigenvectors.
    void
    printSolution ( const Epetra_Vector & x,
                    const std::string   & filenameAppendix
                  ) const;

    void
    setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper );

    // This function is necessary to break the circular dependency with the
    // LOCA_Stepper object to allow for a clean termination
    void
    releaseLocaStepper();

    // needed, e.g., to allocate the linear system in the main file
    virtual Teuchos::RCP<const Epetra_Map>
    getMap() const;

    virtual Teuchos::RCP<Epetra_Vector>
    createSystemVector( const Ginla::State & state ) const;

    Teuchos::RCP<Ginla::State>
    createState(  const Epetra_Vector & x ) const;

private:

    enum continuationType
    {
        ONEPARAMETER,
        TURNINGPOINT
    };
    
    Teuchos::RCP<Epetra_Vector>
    createX(  const Ginla::State & state );

    int
    numDigits ( const int i );

    //! Print method for the continuation in one parameter.
    void
    printSolutionOneParameterContinuation ( const Teuchos::RCP<const Ginla::State> & state );

    //! Print method for turning point continuation continuation.
    void
    printSolutionTurningPointContinuation ( const Teuchos::RCP<const Ginla::State> & state );

    //! Write statistics about the current continuation step to the file
    //! \c outputDataFileName_ .
    void
    writeContinuationStats ( const Teuchos::RCP<const Ginla::State> & state );

    //! Translate an Epetra_Comm into a Teuchos::Comm<int>,
    //! no matter the Thyra::Ordinal.
    Teuchos::RCP<const Teuchos::Comm<int> >
    create_CommInt ( const Teuchos::RCP<const Epetra_Comm> &epetraComm );


private:
    const Teuchos::RCP<Ginla::Operator::Virtual> glOperator_;
    const Teuchos::RCP<Ginla::Perturbation::Virtual> perturbation_;
    Teuchos::RCP<Epetra_CrsMatrix> preconditioner_;

    Teuchos::RCP<const LOCA::Stepper> stepper_;

    const Teuchos::RCP<Ginla::Komplex> komplex_;
    
    const Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter_;
    const Teuchos::RCP<Ginla::IO::StateWriter> stateWriter_;
    
    bool firstTime_;
    
    unsigned int maxNumDigits_; 
    
    continuationType continuationType_;
};

  } // namespace LocaSystem
} // namespace GL

#endif // GINLA_LOCASYSTEM_DEFAULT_H
