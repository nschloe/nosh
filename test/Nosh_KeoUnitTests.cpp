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

#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Vector.h>
#include <LOCA_Parameter_Vector.H>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Constant.hpp"
#include "nosh/VectorField_ExplicitValues.hpp"
#include "nosh/ParameterMatrix_Keo.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace
{

// =============================================================================
void
testKeo(const std::string & inputFileNameBase,
        const double initMu,
        const double controlNormOne,
        const double controlNormInf,
        const double controlSum,
        const double controlSumReal,
        Teuchos::FancyOStream & out,
        bool & success)
{
  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> eComm =
    Teuchos::rcp<Epetra_MpiComm> (new Epetra_MpiComm (MPI_COMM_WORLD));
#else
  Teuchos::RCP<Epetra_SerialComm> eComm =
    Teuchos::rcp<Epetra_SerialComm> (new Epetra_SerialComm());
#endif

  std::string inputFileName = "data/" + inputFileNameBase + ".e";
  // =========================================================================
  // Read the data from the file.
  Teuchos::RCP<Nosh::StkMesh> mesh =
    Teuchos::rcp(new Nosh::StkMesh(eComm, inputFileName, 0));

  // Cast the data into something more accessible.
  Teuchos::RCP<Epetra_Vector> z =
    mesh->createComplexVector("psi");

  std::cout << mesh->getComplexNonOverlapMap() << std::endl;

  Teuchos::RCP<Nosh::VectorField::Virtual> mvp =
    Teuchos::rcp(new Nosh::VectorField::ExplicitValues(*mesh, "A", initMu));

  std::cout << mesh->getComplexNonOverlapMap() << std::endl;

  // Set the thickness field.
  Teuchos::RCP<Nosh::ScalarField::Virtual> thickness =
    Teuchos::rcp(new Nosh::ScalarField::Constant(*mesh, 1.0));

  std::cout << mesh->getComplexNonOverlapMap() << std::endl;

  std::cout << "a" << std::endl;
  Nosh::ParameterMatrix::Keo keo(mesh, thickness, mvp);
  std::cout << "b" << std::endl;

  // Explicitly create the kinetic energy operator.
  std::map<std::string, double> params;
  params["mu"] = initMu;

  keo.refill(params);

  // Compute matrix norms as hashes.
  // Don't check for NormFrobenius() as this one doesn't work for matrices
  // with overlapping maps.
  double normOne = keo.NormOne();
  double normInf = keo.NormInf();

  // check the values
  TEST_FLOATING_EQUALITY(normOne, controlNormOne, 1.0e-12);
  TEST_FLOATING_EQUALITY(normInf, controlNormInf, 1.0e-12);

  const Epetra_Map & map = keo.DomainMap();
  double sum;
  Epetra_Vector u(map);
  Epetra_Vector Ku(map);

  // Add up all the entries of the matrix.
  TEUCHOS_ASSERT_EQUALITY(0, u.PutScalar(1.0));
  TEUCHOS_ASSERT_EQUALITY(0, keo.Apply(u, Ku));
  TEUCHOS_ASSERT_EQUALITY(0, u.Dot(Ku, &sum));
  TEST_FLOATING_EQUALITY(sum, controlSum, 1.0e-10);

  // Sum over all the "real parts" of the matrix.
  // Remember that a 2x2 block corresponding to z is composed as
  // [ Re(z) -Im(z) ]
  // [ Im(z)  Re(z) ].
  // Build vector [ 1, 0, 1, 0, ... ]:
  double one  = 1.0;
  double zero = 0.0;
  for (int k = 0; k < map.NumMyPoints(); k++) {
    if (map.GID(k) % 2 == 0)
      u.ReplaceMyValues(1, &one, &k);
    else
      u.ReplaceMyValues(1, &zero, &k);
  }
  TEUCHOS_ASSERT_EQUALITY(0, keo.Apply(u, Ku));
  TEUCHOS_ASSERT_EQUALITY(0, u.Dot(Ku, &sum));
  TEST_FLOATING_EQUALITY(sum, controlSumReal, 1.0e-10);

  // Sum over all the "imaginary parts" of the matrix.
  // Build vector [ 0, 1, 0, 1, ... ]:
  Epetra_Vector v(map);
  for (int k = 0; k < map.NumMyPoints(); k++) {
    if (map.GID(k) % 2 == 0)
      v.ReplaceMyValues(1, &zero, &k);
    else
      v.ReplaceMyValues(1, &one, &k);
  }
  TEUCHOS_ASSERT_EQUALITY(0, keo.Apply(u, Ku));
  TEUCHOS_ASSERT_EQUALITY(0, v.Dot(Ku, &sum));
  // The matrix is Hermitian, so just test that the sum of
  // the imaginary parts is (close to) 0.
  // Don't use TEST_FLOATING_EQUALITY as this one checks
  // the *relative* error.
  TEST_COMPARE(fabs(sum), <, 1.0e-12);

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
