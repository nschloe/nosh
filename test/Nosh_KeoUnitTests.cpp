// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2011  Nico Schl\"omer
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
#include <map>
#include <string>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Teuchos_ParameterList.hpp>

#include <nosh.hpp>

#include <Teuchos_UnitTestHarness.hpp>

namespace
{

// =============================================================================
void
testKeo(
    const std::string & inputFileNameBase,
    const double initMu,
    const double controlNormOne,
    const double controlNormInf,
    const double controlSum,
    const double controlSumReal,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  std::string inputFileName = "data/" + inputFileNameBase + ".e";

  // Read the data from the file.
  auto mesh = Nosh::read(inputFileName);

  // Cast the data into something more accessible.
  auto z = mesh->getComplexVector("psi");

  auto mvp = std::make_shared<Nosh::VectorField::ExplicitValues>(*mesh, "A", initMu);
  auto thickness = std::make_shared<Nosh::ScalarField::Constant>(*mesh, 1.0);

  Nosh::ParameterMatrix::Keo keo(mesh, thickness, mvp);

  // Explicitly create the kinetic energy operator.
  keo.setParameters({{"mu", initMu}});

  // TODO revive these tests
  //// Compute matrix norms as hashes.
  //// Don't check for NormFrobenius() as this one doesn't work for matrices with
  //// overlapping maps.
  //TEST_FLOATING_EQUALITY(
  //    keo.NormOne(),
  //    controlNormOne,
  //    1.0e-12
  //    );
  //TEST_FLOATING_EQUALITY(
  //    keo.NormInf(),
  //    controlNormInf,
  //    1.0e-12
  //    );

  auto map = keo.getDomainMap();
  Tpetra::Vector<double,int,int> u(map);
  Tpetra::Vector<double,int,int> Ku(map);

  // Add up all the entries of the matrix.
  u.putScalar(1.0);
  keo.apply(u, Ku);
  const double sum = u.dot(Ku);
  TEST_FLOATING_EQUALITY(sum, controlSum, 1.0e-10);

  // Sum over all the "real parts" of the matrix.
  // Remember that a 2x2 block corresponding to z is composed as
  // [ Re(z) -Im(z) ]
  // [ Im(z)  Re(z) ].
  // Build vector [ 1, 0, 1, 0, ... ]:
  double one  = 1.0;
  double zero = 0.0;
  for (size_t k = 0; k < map->getNodeNumElements(); k++) {
    if (map->getGlobalElement(k) % 2 == 0) {
      u.replaceLocalValue(k, one);
    } else {
      u.replaceLocalValue(k, zero);
    }
  }
  keo.apply(u, Ku);
  TEST_FLOATING_EQUALITY(
      u.dot(Ku),
      controlSumReal,
      1.0e-10
      );

  // Sum over all the "imaginary parts" of the matrix.
  // Build vector [ 0, 1, 0, 1, ... ]:
  Tpetra::Vector<double,int,int> v(map);
  for (size_t k = 0; k < map->getNodeNumElements(); k++) {
    if (map->getGlobalElement(k) % 2 == 0) {
      u.replaceLocalValue(k, zero);
    } else {
      u.replaceLocalValue(k, one);
    }
  }
  keo.apply(u, Ku);
  // The matrix is Hermitian, so just test that the sum of the imaginary parts
  // is (close to) 0.
  // Don't use TEST_FLOATING_EQUALITY as this one checks the *relative* error.
  TEST_COMPARE(
      fabs(v.dot(Ku)),
      <,
      1.0e-12
      );

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(Nosh, KeoRectangleSmallHashes)
//{
//  std::string inputFileNameBase = "rectanglesmall";
//
//  double mu = 1.0e-2;
//  double controlNormOne = 10.224658806561596;
//  double controlNormInf = controlNormOne;
//  double controlSumReal = 0.0063121712308067401;
//  double controlSum     = 2 * controlSumReal;
//
//  testKeo(inputFileNameBase,
//          mu,
//          controlNormOne,
//          controlNormInf,
//          controlSum,
//          controlSumReal,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, KeoPacmanHashes)
{
  std::string inputFileNameBase = "pacman";

  double mu = 1.0e-2;
  double controlNormOne = 10.000520856079092;
  double controlNormInf = controlNormOne;
  double controlSumReal = 0.37044264296585938;
  double controlSum     = 2 * controlSumReal;

  testKeo(inputFileNameBase,
          mu,
          controlNormOne,
          controlNormInf,
          controlSum,
          controlSumReal,
          out,
          success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(Nosh, KeoCubeSmallHashes)
//{
//  std::string inputFileNameBase = "cubesmall";
//
//  double mu = 1.0e-2;
//  double controlNormOne = 10.058364522531498;
//  double controlNormInf = controlNormOne;
//  double controlSumReal = 8.3541623155714007e-05;
//  double controlSum     = 2 * controlSumReal;
//
//  testKeo(inputFileNameBase,
//          mu,
//          controlNormOne,
//          controlNormInf,
//          controlSum,
//          controlSumReal,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, KeoBrickWHoleHashes)
{
  std::string inputFileNameBase = "brick-w-hole";

  double mu = 1.0e-2;
  double controlNormOne = 15.131119904340618;
  double controlNormInf = controlNormOne;
  double controlSumReal = 0.16763276012920181;
  double controlSum     = 2 * controlSumReal;

  testKeo(inputFileNameBase,
          mu,
          controlNormOne,
          controlNormInf,
          controlSum,
          controlSumReal,
          out,
          success);
}
// ============================================================================
} // namespace
