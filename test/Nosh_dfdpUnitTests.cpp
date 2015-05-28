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
#include <LOCA_Parameter_Vector.H>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Constant.hpp"
#include "nosh/ParameterMatrix_Keo.hpp"
#include "nosh/ParameterMatrix_DKeoDP.hpp"
#include "nosh/VectorField_ExplicitValues.hpp"
#include "nosh/ModelEvaluator_Nls.hpp"

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

namespace
{
// ===========================================================================
void
computeFiniteDifference_(
  const EpetraExt::ModelEvaluator & modelEval,
  EpetraExt::ModelEvaluator::InArgs & inArgs,
  EpetraExt::ModelEvaluator::OutArgs & outArgs,
  const std::shared_ptr<const Epetra_Vector> & p,
  const int paramIndex,
  const std::shared_ptr<Epetra_Vector> & fdiff
)
{
  const double eps = 1.0e-6;
  Teuchos::RCP<Epetra_Vector> pp = Teuchos::rcp(new Epetra_Vector(*p));

  // Store the original parameter value.
  const double origValue = (*p)[paramIndex];

  // Get vector at x-eps.
  (*pp)[paramIndex] = origValue - eps;
  inArgs.set_p(0, pp);
  Teuchos::RCP<Epetra_Vector> f0 =
    Teuchos::rcp(new Epetra_Vector(fdiff->Map()));
  outArgs.set_f(f0);
  modelEval.evalModel(inArgs, outArgs);

  // Get vector at x+eps.
  (*pp)[paramIndex] = origValue + eps;
  inArgs.set_p(0, pp);
  outArgs.set_f(Teuchos::rcp(fdiff));
  modelEval.evalModel(inArgs, outArgs);

  // Calculate the finite difference approx for df/dp.
  TEUCHOS_ASSERT_EQUALITY(0, fdiff->Update(-1.0, *f0, 1.0));
  TEUCHOS_ASSERT_EQUALITY(0, fdiff->Scale(0.5/eps));

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

  Nosh::ModelEvaluator::Nls modelEval(
      mesh,
      keoBuilder,
      DKeoDPBuilder,
      sp,
      1.0,
      thickness,
      z
      );

  // -------------------------------------------------------------------------
  // Perform the finite difference test for all parameters present in the
  // system.
  // Get a finite-difference approximation of df/dp.
  EpetraExt::ModelEvaluator::InArgs inArgs = modelEval.createInArgs();
  inArgs.set_x(Teuchos::rcp(z));
  EpetraExt::ModelEvaluator::OutArgs outArgs = modelEval.createOutArgs();

  // Get a the initial parameter vector.
  const Teuchos::RCP<const Epetra_Vector> p = modelEval.get_p_init(0);
  std::shared_ptr<Epetra_Vector> fdiff(new Epetra_Vector(z->Map()));
  std::shared_ptr<Epetra_Vector> dfdp(new Epetra_Vector(z->Map()));
  EpetraExt::ModelEvaluator::DerivativeMultiVector deriv;
  std::vector<int> paramIndices(1);
  const std::shared_ptr<Epetra_Vector> nullV = nullptr;
  double r;
  // Only test the first parameter "g" for now since we have to set the DKeoDP
  // above. Alternative: Use "mu" above, and take paramIndex 1 here.
  // TODO test both parameters
  //for (int paramIndex = 0; paramIndex < p->GlobalLength(); paramIndex++) {
  for (int paramIndex = 0; paramIndex < 1; paramIndex++) {
    std::cout << paramIndex << std::endl;
    // Get finite difference.
    computeFiniteDifference_(
        modelEval,
        inArgs,
        outArgs,
        Teuchos::get_shared_ptr(p),
        paramIndex,
        fdiff
        );

    // Get the actual derivative.
    inArgs.set_p(0, p);
    paramIndices[0] = paramIndex;
    deriv = EpetraExt::ModelEvaluator::DerivativeMultiVector(
        Teuchos::rcp(dfdp),
        EpetraExt::ModelEvaluator::DERIV_MV_BY_COL,
        paramIndices
        );
    outArgs.set_DfDp(0, deriv);
    outArgs.set_f(Teuchos::rcp(nullV));
    modelEval.evalModel(inArgs, outArgs);

    // Compare the two.
    TEUCHOS_ASSERT_EQUALITY(0, fdiff->Update(-1.0, *dfdp, 1.0));
    TEUCHOS_ASSERT_EQUALITY(0, fdiff->NormInf(&r));
    TEST_COMPARE(r, <, 1.0e-8);
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
