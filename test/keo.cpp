#include <catch.hpp>

#include <map>
#include <string>

#include <Teuchos_DefaultComm.hpp>

#include <nosh.hpp>

// =============================================================================
void
testKeo(
    const std::string & input_filename_base,
    const double initMu,
    const double control_norm_1,
    const double control_norm_inf,
    const double control_sum,
    const double control_sum_real
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
  keo.set_parameters({{"mu", initMu}}, {});

#if 0
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
#endif

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
  REQUIRE(u.dot(Ku) == Approx(control_sum));

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
  REQUIRE(u.dot(Ku) == Approx(control_sum_real));

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
  REQUIRE(fabs(v.dot(Ku)) == Approx(0.0));

  return;
}
// ===========================================================================
#if 0
TEST_CASE("KEO for rectangle mesh", "[rectangle]")
{
  // For reference: The expected KEO matrix:
  //
  // 5.05         0            -2.00604e-16 0            -4.99844    -0.124987  -0.0499844  0.00124987
  // 0            5.05         0            -2.00604e-16 0.124987    -4.99844   -0.00124987 -0.0499844
  // -2.00604e-16 0            5.05         0            -0.0499844  0.00124987 -4.99844    -0.124987
  // 0            -2.00604e-16 0            5.05         -0.00124987 -0.0499844 0.124987    -4.99844
  // -4.99844     0.124987     -0.0499844   -0.00124987  5.05         0         0           0
  // -0.124987    -4.99844     0.00124987   -0.0499844   0            5.05      0           0
  // -0.0499844   -0.00124987  -4.99844     0.124987     0            0         5.05        0
  // 0.00124987   -0.0499844   -0.124987    -4.99844     0            0         0           5.05
  //
  testKeo(
      "rectanglesmall",
      1.0e-2,
      10.224658806561596,
      10.224658806561596,
      0.01262434246161348, // 2* 0.0063121712308067401
      0.0063121712308067401
      );
}
#endif
// ============================================================================
TEST_CASE("KEO for pacman mesh", "[pacman]")
{
  testKeo(
      "pacman",
      1.0e-2,
      10.000520856079092,
      10.000520856079092,
      0.7408852859317188, // 2 * 0.37044264296585938
      0.37044264296585938
      );
}
// ============================================================================
#if 0
TEST_CASE("KEO for cube mesh", "[cube]")
{
  testKeo(
      "cubesmall",
      1.0e-2,
      10.058364522531498,
      10.058364522531498,
      1.67083246311428e-4, // 2 * 8.3541623155714007e-05
      8.3541623155714007e-05
      );
}
#endif
// ============================================================================
TEST_CASE("KEO for brick mesh", "[brick]")
{
  testKeo(
      "brick-w-hole",
      1.0e-2,
      15.131119904340618,
      15.131119904340618,
      0.3352655202584036,
      0.16763276012920181
      );
}
// ============================================================================
