// @HEADER
//
//    Example function calls.
//    Copyright (C) 2012--2014  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Vector.h>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Constant.hpp"
#include "nosh/ParameterMatrix_Keo.hpp"
#include "nosh/ParameterMatrix_DKeoDP.hpp"
#include "nosh/ParameterMatrix_Laplace.hpp"
#include "nosh/VectorField_ExplicitValues.hpp"
#include "nosh/VectorField_ConstantCurl.hpp"
#include "nosh/ModelEvaluator_Nls.hpp"
#include "nosh/ModelEvaluator_Bordered.hpp"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  // Create a communicator for Epetra objects.
#ifdef HAVE_MPI
  std::shared_ptr<const Epetra_MpiComm> eComm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  std::shared_ptr<const Epetra_SerialComm> eComm(new Epetra_SerialComm());
#endif

  // Create output stream. (Handy for multicore output.)
  const Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();

  bool success = true;
  try {
    // ===========================================================================
    // Handle command line arguments.
    Teuchos::CommandLineProcessor myClp;

    myClp.setDocString(
      "Example usage of the model evaluator.\n"
    );

    std::string dataFile = "";
    myClp.setOption ("mesh-file", &dataFile, "Mesh file", true );

    // Print warning for unrecognized arguments and make sure to throw an
    // exception if something went wrong.
    myClp.recogniseAllOptions ( true );
    myClp.throwExceptions(true);

    // Finally, parse the command line.
    myClp.parse(argc, argv);
    // =======================================================================
    // Read the data from the file.
    const int step = 0;
    std::shared_ptr<Nosh::StkMesh> mesh;
    const Teuchos::RCP<Teuchos::Time> readTime =
      Teuchos::TimeMonitor::getNewTimer("Read mesh");
    {
      Teuchos::TimeMonitor tm(*readTime);
      mesh = std::make_shared<Nosh::StkMesh>(eComm, dataFile, step);
    }

    // Cast the data into something more accessible.
    auto psi = mesh->getComplexVector("psi");

    // Set the thickness field.
    auto thickness = std::make_shared<Nosh::ScalarField::Constant>(*mesh, 1.0);

    // - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -
    // Some alternatives for the positive-definite operator.
    // (a) -\Delta (Laplace operator with Neumann boundary)
    //const std::shared_ptr<Nosh::ParameterMatrix::Virtual> keoBuilder =
    //  rcp(new Nosh::ParameterMatrix::Laplace(mesh, thickness));

    // (b) (-i\nabla-A)^2 (Kinetic energy of a particle in magnetic field)
    // (b1) 'A' explicitly given in file.
    const double initMu = 0.0;
    auto mvp = std::make_shared<Nosh::VectorField::ExplicitValues>(
        *mesh, "A", initMu
        );
    const auto keoBuilder = std::make_shared<Nosh::ParameterMatrix::Keo>(
        mesh, thickness, mvp
        );
    const auto DKeoDPBuilder = std::make_shared<Nosh::ParameterMatrix::DKeoDP>(
        mesh, thickness, mvp, "mu"
        );

    // (b2) 'A' analytically given (here with constant curl).
    //      Optionally add a rotation axis u. This is important
    //      if continuation happens as a rotation of the vector
    //      field around an axis.
    //const std::shared_ptr<DoubleVector> b = rcp(new DoubleVector(3));
    //std::shared_ptr<Teuchos::SerialDenseVector<int,double> > u = Teuchos::null;
    //if ( piroParams->isSublist("Rotation vector") )
    //{
    //    u = rcp(new Teuchos::SerialDenseVector<int,double>(3));
    //    Teuchos::ParameterList & rotationVectorList =
    //        piroParams->sublist( "Rotation vector", false );
    //    (*u)[0] = rotationVectorList.get<double>("x");
    //    (*u)[1] = rotationVectorList.get<double>("y");
    //    (*u)[2] = rotationVectorList.get<double>("z");
    //}
    //std::shared_ptr<Nosh::VectorField::Virtual> mvp =
    //  rcp(new Nosh::VectorField::ConstantCurl(mesh, b, u));
    //const std::shared_ptr<Nosh::ParameterMatrix::Virtual> keoBuilder =
    //  rcp(new Nosh::ParameterMatrix::Keo(mesh, thickness, mvp));
    // (b3) 'A' analytically given in a class you write yourself, derived
    //      from Nosh::ParameterMatrix::Virtual.
    // [...]
    // - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -
    // Setup the scalar potential V.
    // (a) A constant potential.
    auto sp = std::make_shared<Nosh::ScalarField::Constant>(*mesh, -1.0);

    //const double T = 0.0;
    // (b) One you built yourself by deriving from Nosh::ScalarField::Virtual.
    //std::shared_ptr<Nosh::ScalarField::Virtual> sp =
    //rcp(new MyScalarField(mesh));
    // - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

    // Finally, create the model evaluator.
    // This is the most important object in the whole stack.
    const double g = 1.0;
    std::shared_ptr<Nosh::ModelEvaluator::Virtual> nlsModel;
    const Teuchos::RCP<Teuchos::Time> meTime =
      Teuchos::TimeMonitor::getNewTimer("Create model evaluator");
    {
      Teuchos::TimeMonitor tm(*meTime);
      nlsModel = std::make_shared<Nosh::ModelEvaluator::Nls>(
            mesh,
            keoBuilder,
            DKeoDPBuilder,
            sp,
            g,
            thickness,
            psi
            );
    }

    std::shared_ptr<Nosh::ModelEvaluator::Virtual> modelEvaluator;
    const bool useBordering = false;
    if (useBordering) {
      // Use i*psi as bordering.
      std::shared_ptr<Epetra_Vector> bordering(new Epetra_Vector(psi->Map()));
      for (int k = 0; k < psi->Map().NumMyElements()/2; k++) {
        (*bordering)[2*k] = - (*psi)[2*k+1];
        (*bordering)[2*k+1] = (*psi)[2*k];
        //(*bordering)[2*k]   = 1.0;
        //(*bordering)[2*k+1] = 0.0;
      }
      //bordering->Random();
      // Initial value for the extra variable.
      double lambda = 0.0;
      modelEvaluator = std::make_shared<Nosh::ModelEvaluator::Bordered>(
          nlsModel,
          bordering,
          lambda
          );
    } else {
      modelEvaluator = nlsModel;
    }

    const Teuchos::RCP<Teuchos::Time> fxTime =
      Teuchos::TimeMonitor::getNewTimer("F(x)");
    EpetraExt::ModelEvaluator::InArgs inArgs = modelEvaluator->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = modelEvaluator->createOutArgs();
    Teuchos::RCP<Epetra_Operator> nullOp = Teuchos::null;

    // Evaluate the nonlinear Schroedinger equation F(X).
    // Set the parameter vector.
    Teuchos::RCP<const Epetra_Vector> p = modelEvaluator->get_p_init(0);
    inArgs.set_p(0, p);
    // Set the X in inArgs.
    inArgs.set_x(Teuchos::rcp(psi));
    // Set an empty FX in outArgs. This will be filled after the call to evalModel.
    std::shared_ptr<Epetra_Vector> fx(
        new Epetra_Vector(*modelEvaluator->get_f_map())
        );
    outArgs.set_f(Teuchos::rcp(fx));
    modelEvaluator->evalModel(inArgs, outArgs);
    {
      Teuchos::TimeMonitor tm(*fxTime);
      // Fill outArgs.
      modelEvaluator->evalModel(inArgs, outArgs);
      // fx is now filled with F(X).
    }
    //std::cout << *fx << std::endl;
    // Reset to null to make sure it's not refilled the next time evalModel is called.
    Teuchos::RCP<Epetra_Vector> null = Teuchos::null;
    outArgs.set_f(null);

    // Get the Jacobian.
    const Teuchos::RCP<Teuchos::Time> getJTime =
      Teuchos::TimeMonitor::getNewTimer("Get Jacobian");
    Teuchos::RCP<Epetra_Operator> jac;
    {
      Teuchos::TimeMonitor tm(*getJTime);
      jac = modelEvaluator->create_W();
      outArgs.set_W(jac);
      modelEvaluator->evalModel(inArgs, outArgs);
      // Now, jac contains the Jacobian operator.
      outArgs.set_W(nullOp);
    }

    // Apply Jacobian
    Epetra_Vector X(jac->OperatorDomainMap());
    X.Random();
    Epetra_Vector Y(jac->OperatorRangeMap());
    const Teuchos::RCP<Teuchos::Time> applyJTime =
      Teuchos::TimeMonitor::getNewTimer("Apply Jacobian");
    {
      Teuchos::TimeMonitor tm(*applyJTime);
      jac->Apply(X, Y);
    }

    // Get the preconditioner.
    const Teuchos::RCP<Teuchos::Time> getPTime =
      Teuchos::TimeMonitor::getNewTimer("Get Preconditioner");
    Teuchos::RCP<Epetra_Operator> prec;
    {
      Teuchos::TimeMonitor tm(*getPTime);
      prec = modelEvaluator->create_WPrec()->PrecOp;
      outArgs.set_WPrec(prec);
      modelEvaluator->evalModel(inArgs, outArgs);
      // Now, prec contains the Jacobian operator.
      outArgs.set_WPrec(nullOp);
    }

    // Apply preconditioner
    Epetra_Vector X2(prec->OperatorDomainMap());
    X2.Random();
    Epetra_Vector Y2(prec->OperatorRangeMap());
    const Teuchos::RCP<Teuchos::Time> applyPTime =
      Teuchos::TimeMonitor::getNewTimer("Apply Preconditioner");
    {
      Teuchos::TimeMonitor tm(*applyPTime);
      prec->Apply(X2, Y2);
    }

    // Write out data.
    const Teuchos::RCP<Teuchos::Time> writeTime =
      Teuchos::TimeMonitor::getNewTimer("Write");
    {
      Teuchos::TimeMonitor tm(*writeTime);
      mesh->openOutputChannel(".", "output");
      // TODO how to write out psi?
      mesh->write(0.0);
    }

    // Print timing data.
    Teuchos::TimeMonitor::summarize();
  } catch (Teuchos::CommandLineProcessor::HelpPrinted) {
  } catch (Teuchos::CommandLineProcessor::ParseError) {
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
