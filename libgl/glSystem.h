/*! Jacobian system for the Ginzburg--Landau problem.
 *  This routine can be used as an interface to NOX.
 ******************************************************************************/
#ifndef GLSYSTEM_H
#define GLSYSTEM_H

#include "ginzburgLandau.h"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>

#include <Epetra_CrsMatrix.h>

#include <Teuchos_ParameterList.hpp>

#include <NOX_Epetra_Interface_Required.H> // NOX base class
#include <NOX_Epetra_Interface_Jacobian.H> // NOX base class
#include <LOCA_Epetra_Interface_Required.H> // LOCA base class
#include <LOCA_Epetra_Interface_TimeDependent.H> // LOCA base class
#include <NOX_Abstract_PrePostOperator.H>

#include <LOCA_Parameter_Vector.H>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <NOX_Abstract_Group.H>

#include <LOCA_Stepper.H>

typedef Tpetra::Vector<double_complex, int> ComplexVector;

class GlSystem: public NOX::Epetra::Interface::Jacobian,
		public LOCA::Epetra::Interface::TimeDependent {
public:

	//! Constructor with initial guess.
	GlSystem( GinzburgLandau::GinzburgLandau &gl,
			const Teuchos::RCP<const Epetra_Comm> eComm,
			const Teuchos::RCP<ComplexVector> psi,
            const std::string outputDir = "data",
            const std::string outputDataFileName = "continuationData.dat",
            const std::string outputFileFormat = "VTK",
            const std::string solutionFileNameBase = "solutionStep",
            const std::string nullvectorFileNameBase = "nullvectorStep" );

	// Constructor without initial guess.
	GlSystem( GinzburgLandau::GinzburgLandau &gl,
			const Teuchos::RCP<const Epetra_Comm> eComm,
            const std::string outputDir = "data",
            const std::string outputDataFileName = "continuationData.dat",
            const std::string outputFileFormat = "VTK",
            const std::string solutionFileNameBase = "solutionStep",
            const std::string nullvectorFileNameBase = "nullvectorStep" );

	//! Destructor
	~GlSystem();

	//! Evaluate the Ginzburg--Landau functions at a given state defined
	//! by the input vector x.
	virtual bool
	computeF(const Epetra_Vector &x, Epetra_Vector &F,
			const NOX::Epetra::Interface::Required::FillType fillFlag =
					Residual);

	//! Evaluate the Jacobian matrix of the Ginzburg--Landau problem
	//! at a given state defined by the input vector x.
	virtual bool
	computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

	virtual bool
	computeShiftedMatrix(double alpha, double beta, const Epetra_Vector &x,
			Epetra_Operator &A);

	//! Dummy preconditioner function. So far does nothing but throwing
	//! an exception when called.
	virtual bool
	computePreconditioner(const Epetra_Vector &x, Epetra_Operator &Prec,
			Teuchos::ParameterList* precParams = 0) const;

	//! Returns the current state. Not necessarily a solution to the problem.
	//! @return Reference-counted pointer to the current state.
	Teuchos::RCP<Epetra_Vector>
	getSolution() const;

	//! Returns the current Jacobian.
	//! @return Reference-counted pointer to the Jacobian.
	Teuchos::RCP<Epetra_CrsMatrix>
	getJacobian() const;

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

	//! Explicitly print the solution x along with the problem parameters
	//! to the file fileName.
	void
	writeSolutionToFile( const Epetra_Vector &x,
	                     const std::string &filePath) const;

	void
	writeAbstractStateToFile( const Epetra_Vector &x,
	                          const std::string &filePath) const;

	Teuchos::RCP<Epetra_Vector>
	getGlSystemVector( const Teuchos::RCP<const ComplexVector> psi ) const;

private:

	int
	realIndex2complexIndex(const int realIndex) const;

	void
	real2complex(const Epetra_Vector &x, ComplexVector &psi) const;

	void
	complex2real(const ComplexVector &psi, Epetra_Vector &x) const;

	void
	makeRealMap(const Teuchos::RCP<const Tpetra::Map<int> > complexMap);

	enum jacCreator {
		ONLY_GRAPH, VALUES
	};

	enum continuationType {
		ONEPARAMETER,
		TURNINGPOINT
	};

	bool
	createJacobian(const jacCreator jc, const Epetra_Vector &x);

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

	continuationType continuationType_;

	int NumRealUnknowns_;
	int NumMyElements_;
	int NumComplexUnknowns_;

	Teuchos::RCP<const LOCA::Stepper> stepper_;

	GinzburgLandau::GinzburgLandau Gl_;
	const Teuchos::RCP<const Epetra_Comm> EComm_;
	Teuchos::RCP<const Teuchos::Comm<int> > TComm_;
	Teuchos::RCP<Epetra_Map> RealMap_;
	Teuchos::RCP<Epetra_Map> EverywhereMap_;
	Teuchos::RCP<const Tpetra::Map<int> > ComplexMap_;
	Epetra_Vector *rhs_;
	Teuchos::RCP<Epetra_CrsGraph> Graph_;
	Teuchos::RCP<Epetra_CrsMatrix> jacobian_;
	Teuchos::RCP<Epetra_Vector> initialSolution_;

	std::string outputDir_;
	const std::string solutionFileNameBase_;
	const std::string nullvectorFileNameBase_;
	const std::string outputFileFormat_;
	const std::string outputDataFileName_;
};

#endif // GLSYSTEM_H
