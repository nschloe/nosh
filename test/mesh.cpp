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

#include <nosh.hpp>

#include <Teuchos_UnitTestHarness.hpp>

namespace
{

// =============================================================================
void
testMesh(
    const std::string & input_filename_base,
    const unsigned int control_num_nodes,
    const double control_vol_norm_1,
    const double control_vol_norm_2,
    const double control_vol_norm_inf,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  auto comm =  Teuchos::DefaultComm<int>::getComm();
  const int numProc = comm->getSize();

  // Read the data from the file.
  std::string input_filename;
  if (numProc == 1) {
    input_filename = "data/" + input_filename_base + ".e";
  } else {
    input_filename = "data/" + input_filename_base + "-split.par";
  }

  // Read the data from the file.
  auto mesh = nosh::read(input_filename);

  const unsigned int num_nodes = mesh->map()->getGlobalNumElements();
  TEUCHOS_ASSERT_EQUALITY(num_nodes, control_num_nodes);

  const auto control_vols = mesh->control_volumes();
  TEST_FLOATING_EQUALITY(
      control_vols->norm1(),
      control_vol_norm_1,
      1.0e-12
      );
  TEST_FLOATING_EQUALITY(
      control_vols->norm2(),
      control_vol_norm_2,
      1.0e-12
      );
  TEST_FLOATING_EQUALITY(
      control_vols->normInf(),
      control_vol_norm_inf,
      1.0e-12
      );

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(nosh, MeshRectangleSmallHashes)
//{
//  std::string input_filename_base = "rectanglesmall";
//
//  unsigned int num_nodes = 4;
//  double control_vol_norm_1 = 10.0;
//  double control_vol_norm_2 = 5.0;
//  double control_vol_norm_inf = 2.5;
//
//  testMesh(input_filename_base,
//           num_nodes,
//           control_vol_norm_1,
//           control_vol_norm_2,
//           control_vol_norm_inf,
//           out,
//           success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, MeshPacmanHashes)
{
  std::string input_filename_base = "pacman";

  unsigned int num_nodes = 409;
  double control_vol_norm_1 = 302.5227007210103;
  double control_vol_norm_2 = 15.38575790933914;
  double control_vol_norm_inf = 1.127797467043659;

  testMesh(
      input_filename_base,
      num_nodes,
      control_vol_norm_1,
      control_vol_norm_2,
      control_vol_norm_inf,
      out,
      success
      );
}
// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, MeshShellHashes)
//{
//  std::string input_filename_base = "shell";
//
//  unsigned int num_nodes = 5;
//  double control_vol_norm_1 = 3.46410161513775;
//  double control_vol_norm_2 = 1.63299316185545;
//  double control_vol_norm_inf = 1.15470053837925;
//
//  testMesh(input_filename_base,
//           num_nodes,
//           control_vol_norm_1,
//           control_vol_norm_2,
//           control_vol_norm_inf,
//           out,
//           success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, MeshSphereHashes)
{
  std::string input_filename_base = "sphere";

  unsigned int num_nodes = 82;
  double control_vol_norm_1 = 11.9741927059035;
  double control_vol_norm_2 = 1.39047542328083;
  double control_vol_norm_inf = 0.198927169088121;

  testMesh(
      input_filename_base,
      num_nodes,
      control_vol_norm_1,
      control_vol_norm_2,
      control_vol_norm_inf,
      out,
      success
      );
}
// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, MeshCubeSmallHashes)
//{
//  std::string input_filename_base = "cubesmall";
//
//  unsigned int num_nodes = 8;
//  double control_vol_norm_1 = 10.0;
//  double control_vol_norm_2 = 3.535533905932738;
//  double control_vol_norm_inf = 1.25;
//
//  testMesh(input_filename_base,
//           num_nodes,
//           control_vol_norm_1,
//           control_vol_norm_2,
//           control_vol_norm_inf,
//           out,
//           success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, MeshBrickWHoleHashes)
{
  std::string input_filename_base = "brick-w-hole";

  unsigned int num_nodes = 744;
  double control_vol_norm_1 = 388.686291694641;
  double control_vol_norm_2 = 16.6614019419857;
  double control_vol_norm_inf = 1.46847345474977;

  testMesh(
      input_filename_base,
      num_nodes,
      control_vol_norm_1,
      control_vol_norm_2,
      control_vol_norm_inf,
      out,
      success
      );
}
// ============================================================================
} // namespace
