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
    const std::string & input_filename_base,
    const double initMu,
    const double control_norm_1,
    const double control_norm_inf,
    const double control_sum,
    const double control_sum_real,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  // Read the data from the file.
  auto comm =  Teuchos::DefaultComm<int>::getComm();
  const int size = comm->getSize();
  const std::string input_filename = (size == 1) ?
    "data/" + input_filename_base + ".h5m" :
    "data/" + input_filename_base + "-" + std::to_string(size) + ".h5m"
    ;

  // Read the data from the file.
  auto mesh = nosh::read(input_filename);

  // Cast the data into something more accessible.
  auto z = mesh->get_complex_vector("psi");

  auto mvp = std::make_shared<nosh::vector_field::explicit_values>(*mesh, "A", initMu);
  auto thickness = std::make_shared<nosh::scalar_field::constant>(*mesh, 1.0);

  nosh::parameter_matrix::keo keo(mesh, thickness, mvp);

  // Explicitly create the kinetic energy operator.
  keo.set_parameters({{"mu", initMu}});

  // show me the matrix
  for (int row=0; row < 8; row++) {
    Teuchos::ArrayView<const int> cols;
    Teuchos::ArrayView<const double> vals;
    keo.getLocalRowView(row, cols, vals);
    std::cout << "row: " << row << std::endl;
    for (int k = 0; k < cols.size(); k++) {
      std::cout << "  col: "  << cols[k] << "   val: " << vals[k] << std::endl;
    }
  }

  //// TODO revive these tests
  //// Compute matrix norms as hashes.
  //// Don't check for NormFrobenius() as this one doesn't work for matrices with
  //// overlapping maps.
  (void) control_norm_1;
  (void) control_norm_inf;
  //TEST_FLOATING_EQUALITY(
  //    keo.NormOne(),
  //    control_norm_1,
  //    1.0e-12
  //    );
  //TEST_FLOATING_EQUALITY(
  //    keo.NormInf(),
  //    control_norm_inf,
  //    1.0e-12
  //    );

  auto map = keo.getDomainMap();
  Tpetra::Vector<double,int,int> u(map);
  Tpetra::Vector<double,int,int> Ku(map);

  // Add up all the entries of the matrix.
  u.putScalar(1.0);
  keo.apply(u, Ku);
  TEST_FLOATING_EQUALITY(u.dot(Ku), control_sum, 1.0e-10);

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
      control_sum_real,
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
//TEUCHOS_UNIT_TEST(nosh, KeoRectangleSmallHashes)
//{
//  std::string input_filename_base = "rectanglesmall";
//
//  double mu = 1.0e-2;
//  double control_norm_1 = 10.224658806561596;
//  double control_norm_inf = control_norm_1;
//  double control_sum_real = 0.0063121712308067401;
//  double control_sum     = 2 * control_sum_real;
//  // For reference: The expected KEO matrix:
//  //
//  // 5.05         0            -2.00604e-16 0            -4.99844    -0.124987  -0.0499844  0.00124987
//  // 0            5.05         0            -2.00604e-16 0.124987    -4.99844   -0.00124987 -0.0499844
//  // -2.00604e-16 0            5.05         0            -0.0499844  0.00124987 -4.99844    -0.124987
//  // 0            -2.00604e-16 0            5.05         -0.00124987 -0.0499844 0.124987    -4.99844
//  // -4.99844     0.124987     -0.0499844   -0.00124987  5.05         0         0           0
//  // -0.124987    -4.99844     0.00124987   -0.0499844   0            5.05      0           0
//  // -0.0499844   -0.00124987  -4.99844     0.124987     0            0         5.05        0
//  // 0.00124987   -0.0499844   -0.124987    -4.99844     0            0         0           5.05
//  //
//  testKeo(
//      input_filename_base,
//      mu,
//      control_norm_1,
//      control_norm_inf,
//      control_sum,
//      control_sum_real,
//      out,
//      success
//      );
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, KeoPacmanHashes)
{
  std::string input_filename_base = "pacman";

  double mu = 1.0e-2;
  double control_norm_1 = 10.000520856079092;
  double control_norm_inf = control_norm_1;
  double control_sum_real = 0.37044264296585938;
  double control_sum     = 2 * control_sum_real;

  testKeo(
      input_filename_base,
      mu,
      control_norm_1,
      control_norm_inf,
      control_sum,
      control_sum_real,
      out,
      success
      );
}
// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, KeoCubeSmallHashes)
//{
//  std::string input_filename_base = "cubesmall";
//
//  double mu = 1.0e-2;
//  double control_norm_1 = 10.058364522531498;
//  double control_norm_inf = control_norm_1;
//  double control_sum_real = 8.3541623155714007e-05;
//  double control_sum     = 2 * control_sum_real;
//
//  testKeo(input_filename_base,
//          mu,
//          control_norm_1,
//          control_norm_inf,
//          control_sum,
//          control_sum_real,
//          out,
//          success);
//}
//// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, KeoBrickWHoleHashes)
//{
//  std::string input_filename_base = "brick-w-hole";
//
//  double mu = 1.0e-2;
//  double control_norm_1 = 15.131119904340618;
//  double control_norm_inf = control_norm_1;
//  double control_sum_real = 0.16763276012920181;
//  double control_sum     = 2 * control_sum_real;
//
//  testKeo(input_filename_base,
//          mu,
//          control_norm_1,
//          control_norm_inf,
//          control_sum,
//          control_sum_real,
//          out,
//          success);
//}
//// ============================================================================
} // namespace
