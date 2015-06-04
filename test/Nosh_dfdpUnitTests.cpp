// @HEADER
//
//    Unit test for dF/dp.
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

#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Thyra_EpetraModelEvaluator.hpp>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Constant.hpp"
#include "nosh/ParameterMatrix_Keo.hpp"
#include "nosh/ParameterMatrix_DKeoDP.hpp"
#include "nosh/VectorField_ExplicitValues.hpp"
#include "nosh/ModelEvaluator_Nls.hpp"
#include "nosh/ModelEvaluatorT_Nls.hpp"

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

namespace
{
// ===========================================================================
void
computeFiniteDifference_(
  const Thyra::ModelEvaluator<double> & modelEval,
  const Teuchos::RCP<Thyra::VectorBase<double>> & x,
  Teuchos::RCP<Thyra::VectorBase<double> > & p,
  const int paramIndex,
  const Teuchos::RCP<Thyra::VectorBase<double> > & fdiff
)
{
  const double eps = 1.0e-8;
  Teuchos::RCP<Thyra::VectorBase<double> > pp = p->clone_v();

  const double origValue = Thyra::get_ele(*pp, paramIndex);

  Thyra::ModelEvaluatorBase::InArgs<double> inArgs =
    modelEval.createInArgs();
  inArgs.set_x(x);

  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs =
    modelEval.createOutArgs();

  // Get vector at x-eps.
  Thyra::set_ele(paramIndex, origValue - eps, pp());
  inArgs.set_p(0, pp);
  Teuchos::RCP<Thyra::VectorBase<double> > f0 =
    Thyra::createMember(fdiff->space());
  outArgs.set_f(f0);
  modelEval.evalModel(inArgs, outArgs);

  // Get vector at x+eps.
  Thyra::set_ele(paramIndex, origValue + eps, pp());
  inArgs.set_p(0, pp);
  outArgs.set_f(fdiff);
  modelEval.evalModel(inArgs, outArgs);

  // Calculate the finite difference approx for df/dp.
  Thyra::Vp_StV(fdiff(), -1.0, *f0);
  Thyra::scale(0.5/eps, fdiff());

  return;
}
// =============================================================================
void
testDfdp(const std::string & inputFileNameBase,
         const double mu,
         Teuchos::FancyOStream & out,
         bool & success
        )
{
  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  std::shared_ptr<Epetra_MpiComm> eComm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  std::shared_ptr<Epetra_SerialComm> eComm(new Epetra_SerialComm());
#endif

  std::string inputFileName = "data/" + inputFileNameBase + ".e";
  // =========================================================================
  // Read the data from the file.
  std::shared_ptr<Nosh::StkMesh> mesh(
      new Nosh::StkMesh(eComm, inputFileName, 0)
      );

  // Cast the data into something more accessible.
  std::shared_ptr<Epetra_Vector> z = mesh->createComplexVector("psi");

  // Set the thickness field.
  std::shared_ptr<Nosh::ScalarField::Virtual> thickness(
      new Nosh::ScalarField::Constant(*mesh, 1.0)
      );

  std::shared_ptr<Nosh::VectorField::Virtual> mvp(
      new Nosh::VectorField::ExplicitValues(*mesh, "A", mu)
      );
  const std::shared_ptr<Nosh::ParameterMatrix::Virtual> keoBuilder(
      new Nosh::ParameterMatrix::Keo(mesh, thickness, mvp)
      );
  const std::shared_ptr<Nosh::ParameterMatrix::Virtual> DKeoDPBuilder(
      new Nosh::ParameterMatrix::DKeoDP(mesh, thickness, mvp, "g")
      );

  std::shared_ptr<Nosh::ScalarField::Virtual> sp(
      new Nosh::ScalarField::Constant(*mesh, -1.0)
      );

  //Teuchos::RCP<Nosh::ModelEvaluator::Nls> modelEvalE =
  //  Teuchos::rcp(new Nosh::ModelEvaluator::Nls(
  //        mesh,
  //        keoBuilder,
  //        DKeoDPBuilder,
  //        sp,
  //        1.0,
  //        thickness,
  //        z
  //        ));
  //Teuchos::RCP<Thyra::ModelEvaluator<double> > modelEval =
  //  Thyra::epetraModelEvaluator(modelEvalE, Teuchos::null);

  Teuchos::RCP<Thyra::ModelEvaluator<double>> modelEval =
    Teuchos::rcp(new Nosh::ModelEvaluatorT::Nls(
          mesh,
          keoBuilder,
          DKeoDPBuilder,
          sp,
          1.0,
          thickness,
          z
          ));

  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vectorSpaceX =
    modelEval->get_x_space();
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vectorSpaceF =
    modelEval->get_f_space();

  // -------------------------------------------------------------------------
  // Perform the finite difference test for all parameters present in the
  // system.
  // Get a finite-difference approximation of df/dp.
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs =
    modelEval->createInArgs();
  Teuchos::RCP<Thyra::VectorBase<double>> zT =
    Thyra::create_Vector(Teuchos::rcp(z), vectorSpaceX);
  inArgs.set_x(zT);

  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs =
    modelEval->createOutArgs();

  // create parameter vector
  Teuchos::RCP<Thyra::VectorBase<double> > p = Thyra::createMember(
      modelEval->get_p_space(0)
      );
  Teuchos::RCP<Thyra::VectorBase<double> > fdiff = Thyra::createMember(
      modelEval->get_f_space()
      );

  // Get the actual derivatives.
  inArgs.set_p(0, p);
  Teuchos::RCP<Thyra::MultiVectorBase<double> > dfdp =
    Thyra::createMembers(modelEval->get_f_space(), 2);
  Thyra::ModelEvaluatorBase::Derivative<double> deriv(
      dfdp,
      Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL
      );
  outArgs.set_DfDp(0, deriv);
  modelEval->evalModel(inArgs, outArgs);

  // Only test the first parameter "g" for now since we have to set the DKeoDP
  // above. Alternative: Use "mu" above, and take paramIndex 1 here.
  // TODO test both parameters
  //for (int paramIndex = 0; paramIndex < p->GlobalLength(); paramIndex++) {
  std::vector<int> paramIndices(1);
  double r;
  for (int paramIndex = 0; paramIndex < 1; paramIndex++) {
    // Get finite difference.
    computeFiniteDifference_(
        *modelEval,
        zT,
        p,
        paramIndex,
        fdiff
        );

    // Compare the two.
    Thyra::Vp_StV(fdiff(), -1.0, *dfdp->col(0));
    r = Thyra::norm_inf(*fdiff);
    TEST_COMPARE(r, <, 1.0e-7);
  }
  // -------------------------------------------------------------------------

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(Nosh, DfdpRectangleSmallHashes)
//{
//  const std::string inputFileNameBase = "rectanglesmall";
//  const double mu = 1.0e-2;
//  testDfdp(inputFileNameBase, mu, out, success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, DfdpPacmanHashes)
{
  const std::string inputFileNameBase = "pacman";
  const double mu = 1.0e-2;
  testDfdp(inputFileNameBase, mu, out, success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(Nosh, DfdpCubeSmallHashes)
//{
//  const std::string inputFileNameBase = "cubesmall";
//  const double mu = 1.0e-2;
//  testDfdp(inputFileNameBase, mu, out, success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, DfdpBrickWithHoleHashes)
{
  const std::string inputFileNameBase = "brick-w-hole";
  const double mu = 1.0e-2;
  testDfdp(inputFileNameBase, mu, out, success);
}
// ============================================================================
} // namespace
