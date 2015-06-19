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
#include <string>

#include <Teuchos_ParameterList.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Tpetra_Vector.hpp>

#include <nosh.hpp>

#include <Teuchos_UnitTestHarness.hpp>

namespace
{

// =============================================================================
void
testKeo(
    const std::string & inputFileNameBase,
    const double psiControlNormOne,
    const double psiControlNormInf,
    const std::vector<double> & mvpControlNormsInf,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  //std::string inputFileName = "data/" + inputFileNameBase + ".xmf";
  //Nosh::readXdmf(inputFileName);

  // Read the data from the file.
  std::string inputFileName = "data/" + inputFileNameBase + ".e";
  auto mesh = Nosh::read(inputFileName);

  //return;

  // Cast the data into something more accessible.
  const auto psi = mesh->getComplexVector("psi");
  const auto mvpValues = mesh->getMultiVector("A");

  // Check psi.
  TEST_FLOATING_EQUALITY(
      psi->norm1(),
      psiControlNormOne,
      1.0e-12
      );
  TEST_FLOATING_EQUALITY(
      psi->normInf(),
      psiControlNormInf,
      1.0e-12
      );

  // Check MVP.
  // Only check the infinity-norm here as all other norms only apply to vectors
  // with non-overlapping maps.
  std::vector<double> v(mvpValues->getNumVectors());
  mvpValues->normInf(Teuchos::ArrayView<double>(v));
  TEST_COMPARE_FLOATING_ARRAYS(
      v,
      mvpControlNormsInf,
      1.0e-12
      );

  Teuchos::TimeMonitor::summarize();

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(Nosh, KeoRectangleSmallHashes)
//{
//  std::string inputFileNameBase = "rectanglesmall";
//
//  const double psiControlNormOne = 4.0;
//  const double psiControlNormInf = 1.0;
//  std::vector<double> mvpControlNormsInf(3);
//  mvpControlNormsInf[0] = 0.25;
//  mvpControlNormsInf[1] = 2.5;
//  mvpControlNormsInf[2] = 0.0;
//
//  testKeo(inputFileNameBase,
//          psiControlNormOne,
//          psiControlNormInf,
//          mvpControlNormsInf,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, KeoPacmanHashes)
{
  std::string inputFileNameBase = "pacman";

  const double psiControlNormOne = 409.0;
  const double psiControlNormInf = 1.0;
  std::vector<double> mvpControlNormsInf(3);
  mvpControlNormsInf[0] = 4.999111652374270;
  mvpControlNormsInf[1] = 5.0;
  mvpControlNormsInf[2] = 0.0;

  testKeo(inputFileNameBase,
          psiControlNormOne,
          psiControlNormInf,
          mvpControlNormsInf,
          out,
          success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(Nosh, KeoShellHashes)
//{
//  std::string inputFileNameBase = "shell";
//
//  const double psiControlNormOne = 5.0;
//  const double psiControlNormInf = 1.0;
//  std::vector<double> mvpControlNormsInf(3);
//  mvpControlNormsInf[0] = 0.5;
//  mvpControlNormsInf[1] = 0.5;
//  mvpControlNormsInf[2] = 0.0;
//
//  testKeo(inputFileNameBase,
//          psiControlNormOne,
//          psiControlNormInf,
//          mvpControlNormsInf,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, KeoSphereHashes)
{
  std::string inputFileNameBase = "sphere";

  const double psiControlNormOne = 82.0;
  const double psiControlNormInf = 1.0;
  std::vector<double> mvpControlNormsInf(3);
  mvpControlNormsInf[0] = 0.492403864860535;
  mvpControlNormsInf[1] = 0.468303918838501;
  mvpControlNormsInf[2] = 0.0;

  testKeo(inputFileNameBase,
          psiControlNormOne,
          psiControlNormInf,
          mvpControlNormsInf,
          out,
          success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(Nosh, KeoCubeSmallHashes)
//{
//  std::string inputFileNameBase = "cubesmall";
//
//  const double psiControlNormOne = 8.0;
//  const double psiControlNormInf = 1.0;
//  std::vector<double> mvpControlNormsInf(3);
//  mvpControlNormsInf[0] = 0.25;
//  mvpControlNormsInf[1] = 0.25;
//  mvpControlNormsInf[2] = 0.0;
//
//  testKeo(inputFileNameBase,
//          psiControlNormOne,
//          psiControlNormInf,
//          mvpControlNormsInf,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, KeoBrickWHoleHashes)
{
  std::string inputFileNameBase = "brick-w-hole";

  const double psiControlNormOne = 744.0;
  const double psiControlNormInf = 1.0;
  std::vector<double> mvpControlNormsInf(3);
  mvpControlNormsInf[0] = 2.5;
  mvpControlNormsInf[1] = 2.5;
  mvpControlNormsInf[2] = 0.0;

  testKeo(inputFileNameBase,
          psiControlNormOne,
          psiControlNormInf,
          mvpControlNormsInf,
          out,
          success);
}
// ============================================================================
} // namespace
