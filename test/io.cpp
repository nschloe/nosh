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
    const std::string & input_filename_base,
    const double psi_control_norm_1,
    const double psi_control_norm_inf,
    const std::vector<double> & mvp_control_norms_inf,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  // Read the data from the file.
  auto comm =  Teuchos::DefaultComm<int>::getComm();
  const std::string input_filename = (comm->getSize() == 1) ?
    "data/" + input_filename_base + ".e" :
    "data/" + input_filename_base + "-split.par";

  auto mesh = nosh::read(input_filename);

  //return;

  // Cast the data into something more accessible.
  const auto psi = mesh->get_complex_vector("psi");
  const auto mvpValues = mesh->get_multi_vector("A");

  // Check psi.
  TEST_FLOATING_EQUALITY(
      psi->norm1(),
      psi_control_norm_1,
      1.0e-12
      );
  TEST_FLOATING_EQUALITY(
      psi->normInf(),
      psi_control_norm_inf,
      1.0e-12
      );

  // Check MVP.
  // Only check the infinity-norm here as all other norms only apply to vectors
  // with non-overlapping maps.
  std::vector<double> v(mvpValues->getNumVectors());
  mvpValues->normInf(Teuchos::ArrayView<double>(v));
  TEST_COMPARE_FLOATING_ARRAYS(
      v,
      mvp_control_norms_inf,
      1.0e-12
      );

  Teuchos::TimeMonitor::summarize();

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(nosh, KeoRectangleSmallHashes)
//{
//  std::string input_filename_base = "rectanglesmall";
//
//  const double psi_control_norm_1 = 4.0;
//  const double psi_control_norm_inf = 1.0;
//  std::vector<double> mvp_control_norms_inf(3);
//  mvp_control_norms_inf[0] = 0.25;
//  mvp_control_norms_inf[1] = 2.5;
//  mvp_control_norms_inf[2] = 0.0;
//
//  testKeo(input_filename_base,
//          psi_control_norm_1,
//          psi_control_norm_inf,
//          mvp_control_norms_inf,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, KeoPacmanHashes)
{
  std::string input_filename_base = "pacman";

  const double psi_control_norm_1 = 409.0;
  const double psi_control_norm_inf = 1.0;
  std::vector<double> mvp_control_norms_inf(3);
  mvp_control_norms_inf[0] = 4.999111652374270;
  mvp_control_norms_inf[1] = 5.0;
  mvp_control_norms_inf[2] = 0.0;

  testKeo(input_filename_base,
          psi_control_norm_1,
          psi_control_norm_inf,
          mvp_control_norms_inf,
          out,
          success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, KeoShellHashes)
//{
//  std::string input_filename_base = "shell";
//
//  const double psi_control_norm_1 = 5.0;
//  const double psi_control_norm_inf = 1.0;
//  std::vector<double> mvp_control_norms_inf(3);
//  mvp_control_norms_inf[0] = 0.5;
//  mvp_control_norms_inf[1] = 0.5;
//  mvp_control_norms_inf[2] = 0.0;
//
//  testKeo(input_filename_base,
//          psi_control_norm_1,
//          psi_control_norm_inf,
//          mvp_control_norms_inf,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, KeoSphereHashes)
{
  std::string input_filename_base = "sphere";

  const double psi_control_norm_1 = 82.0;
  const double psi_control_norm_inf = 1.0;
  std::vector<double> mvp_control_norms_inf(3);
  mvp_control_norms_inf[0] = 0.492403864860535;
  mvp_control_norms_inf[1] = 0.468303918838501;
  mvp_control_norms_inf[2] = 0.0;

  testKeo(input_filename_base,
          psi_control_norm_1,
          psi_control_norm_inf,
          mvp_control_norms_inf,
          out,
          success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, KeoCubeSmallHashes)
//{
//  std::string input_filename_base = "cubesmall";
//
//  const double psi_control_norm_1 = 8.0;
//  const double psi_control_norm_inf = 1.0;
//  std::vector<double> mvp_control_norms_inf(3);
//  mvp_control_norms_inf[0] = 0.25;
//  mvp_control_norms_inf[1] = 0.25;
//  mvp_control_norms_inf[2] = 0.0;
//
//  testKeo(input_filename_base,
//          psi_control_norm_1,
//          psi_control_norm_inf,
//          mvp_control_norms_inf,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, KeoBrickWHoleHashes)
{
  std::string input_filename_base = "brick-w-hole";

  const double psi_control_norm_1 = 744.0;
  const double psi_control_norm_inf = 1.0;
  std::vector<double> mvp_control_norms_inf(3);
  mvp_control_norms_inf[0] = 2.5;
  mvp_control_norms_inf[1] = 2.5;
  mvp_control_norms_inf[2] = 0.0;

  testKeo(input_filename_base,
          psi_control_norm_1,
          psi_control_norm_inf,
          mvp_control_norms_inf,
          out,
          success);
}
// ============================================================================
} // namespace
