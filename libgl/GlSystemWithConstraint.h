/*
 * GlSystemWithConstraint.h
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl�mer
 */
#ifndef GLSYSTEMWITHCONSTRAINT_H_
#define GLSYSTEMWITHCONSTRAINT_H_

#include "ginzburgLandau.h"
#include "AbstractStateWriter.h"
#include "GlKomplex.h"
#include "glSystem.h"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

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

class GlSystemWithConstraint: public AbstractStateWriter,
                              public NOX::Epetra::Interface::Jacobian,
                              public NOX::Epetra::Interface::Preconditioner,
                              public LOCA::Epetra::Interface::TimeDependent
{
public:

  //! Constructor with initial guess.
  GlSystemWithConstraint( GinzburgLandau::GinzburgLandau &gl,
                          const Teuchos::RCP<const Epetra_Comm> eComm,
                          const Teuchos::RCP<ComplexVector> psi,
                          const std::string outputDir = "data",
                          const std::string outputDataFileName = "continuationData.dat",
                          const std::string outputFileFormat = "VTK",
                          const std::string solutionFileNameBase = "solutionStep",
                          const std::string nullvectorFileNameBase = "nullvectorStep",
                          const unsigned int maxStepNumberDecimals = 3 );

  // Constructor without initial guess.
  GlSystemWithConstraint( GinzburgLandau::GinzburgLandau &gl,
                          const Teuchos::RCP<const Epetra_Comm> eComm,
                          const std::string outputDir = "data",
                          const std::string outputDataFileName = "continuationData.dat",
                          const std::string outputFileFormat = "VTK",
                          const std::string solutionFileNameBase = "solutionStep",
                          const std::string nullvectorFileNameBase = "nullvectorStep",
                          const unsigned int maxStepNumberDecimals = 3 );

        //! Destructor
        ~GlSystemWithConstraint();

        //! Evaluate the Ginzburg--Landau functions at a given state defined
        //! by the input vector x.
        virtual bool
        computeF(const Epetra_Vector &x, Epetra_Vector &F,
                 const NOX::Epetra::Interface::Required::FillType fillFlag = Residual);

        //! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
        //! at a given state defined by the input vector x.
        virtual bool
        computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

        virtual bool
        computeShiftedMatrix(double alpha,
                             double beta,
                             const Epetra_Vector &x,
                             Epetra_Operator &A );

        //! Dummy preconditioner function. So far does nothing but throwing
        //! an exception when called.
        virtual bool
        computePreconditioner(const Epetra_Vector &x, Epetra_Operator &Prec,
                        Teuchos::ParameterList* precParams = 0);

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
        setParameters(const LOCA::ParameterVector &p);

        //! Set directory to where all output gets written.
        //! @param directory Name of the directory.
        void
        setOutputDir(const string &directory);

        //! Print the solution x along with the continuation parameter conParam
        //! to a file. This function is called internally by LOCA to print
        //! solutions of each continuation step.
        virtual void
        printSolution(const Epetra_Vector &x, double conParam);

        void
        setLocaStepper( const Teuchos::RCP<const LOCA::Stepper> stepper );

        // This function is necessary to break the circular dependency with the
        // LOCA_Stepper object to allow for a clean termination
        void
        releaseLocaStepper();

        //! Explicitly print the solution x along with the problem parameters
        //! to the file fileName.
        void
        writeSolutionToFile( const Epetra_Vector & x,
                             const std::string   & filePath) const;

        void
        writeAbstractStateToFile( const Epetra_Vector & x,
                                  const std::string   & filePath) const;

        // TODO delete?
        const Teuchos::RCP<const GlKomplex>
        getGlKomplex() const;

        // TODO delete
        double
        getH0() const;
        void
        setH0( const double h0 );
        void
        setScaling( const double h0 );
        void
        setChi( const double h0 );


        const Teuchos::RCP<const Epetra_Map>
        getMap() const;

        enum jacCreator {
                ONLY_GRAPH, VALUES
        };

private:

        void
        initialize(const Teuchos::RCP<ComplexVector> psi);

        void
        fillBorderedMatrix( const Teuchos::RCP<      Epetra_CrsMatrix> & extendedMatrix,
                            const Teuchos::RCP<const Epetra_CrsMatrix> & regularMatrix,
                            const Teuchos::Array<double>               & rightBorder,
                                  Teuchos::Array<double>               & lowerBorder,
                                  double                                 d,
                                  bool                                   firstTime
                          ) const;

        int
        PutRow( const Teuchos::RCP<Epetra_CrsMatrix> A,
        		int      Row,
                const int      numIndices,
        		double * values,
        		int    * indices,
        		const bool     firstTime ) const;


        //! Creates a map identical to the input argument \c realMap, but
        //! extended by one entry.
        //! This extra slot is usually used for the phase condition.
        Teuchos::RCP<Epetra_Map>
        createExtendedRealMap( const Epetra_BlockMap & realMap  ) const;

        enum continuationType {
                ONEPARAMETER,
                TURNINGPOINT
        };

        void
        createJacobian( const Epetra_Vector &x );

        //! Print method for the continuation in one parameter.
        void
        printSolutionOneParameterContinuation( const Teuchos::RCP<const ComplexVector> & psi
                                             ) const;

        //! Print method for turning point continuation continuation.
        void
        printSolutionTurningPointContinuation( const Teuchos::RCP<const ComplexVector> & psi
                                         ) const;


        //! Write statistics about the current continuation step to the file
        //! \c outputDataFileName_ .
        void
        writeContinuationStats( const int conStep,
                                const Teuchos::RCP<const ComplexVector> psi ) const;

        //! Translate an Epetra_Comm into a Teuchos::Comm<int>, no matter the Thyra::Ordinal.
        Teuchos::RCP<const Teuchos::Comm<int> >
        create_CommInt( const Teuchos::RCP<const Epetra_Comm> &epetraComm );

        continuationType continuationType_;

        GlSystem::GlSystem glSystem_;

        int NumMyElements_;
        int NumComplexUnknowns_;

        GinzburgLandau::GinzburgLandau Gl_;
        const Teuchos::RCP<const Epetra_Comm> EComm_;
        Teuchos::RCP<const Epetra_BlockMap> regularRealMap_;
        Teuchos::RCP<Epetra_Map> extendedRealMap_;
        Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > ComplexMap_;
        Epetra_Vector *rhs_;
        Teuchos::RCP<Epetra_CrsGraph> Graph_;
        Teuchos::RCP<Epetra_CrsMatrix> jacobian_;
        Teuchos::RCP<Epetra_CrsMatrix> preconditioner_;
        Teuchos::RCP<Epetra_Vector> initialSolution_;

        const std::string solutionFileNameBase_;
        const std::string nullvectorFileNameBase_;
        const std::string outputFileFormat_;
        const std::string outputDataFileName_;

        Teuchos::RCP<GlKomplex> glKomplex_;

        const unsigned int maxStepNumberDecimals_;
        std::string stepNumFileNameFormat_;

        bool firstTime_;
};

#endif /* GLSYSTEMWITHCONSTRAINT_H_ */
