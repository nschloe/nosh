#include <catch.hpp>

#include <string>

#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_Vector.hpp>

#include <nosh.hpp>

// =============================================================================
void
testKeo(
    const std::string & input_filename_base,
    const double psi_control_norm_1,
    const double psi_control_norm_inf,
    const std::vector<double> & mvp_control_norms_inf
    )
{
  // Read the data from the file.
  const auto comm = Teuchos::DefaultComm<int>::getComm();
  const int size = comm->getSize();
  const std::string input_filename = (size == 1) ?
    "data/" + input_filename_base + ".h5m" :
    "data/" + input_filename_base + "-" + std::to_string(size) + ".h5m"
    ;

  auto mesh = nosh::read(input_filename);

  // Cast the data into something more accessible.
  const auto psi = mesh->get_complex_vector("psi");
  const auto mvpValues = mesh->get_multi_vector("A");

  // Check psi.
  REQUIRE(psi->norm1() == Approx(psi_control_norm_1));
  REQUIRE(psi->normInf() == Approx(psi_control_norm_inf));

  // Check MVP.
  // Only check the infinity-norm here as all other norms only apply to vectors
  // with non-overlapping maps.
  std::vector<double> v(mvpValues->getNumVectors());
  mvpValues->normInf(Teuchos::ArrayView<double>(v));
  for (size_t k = 0; k < v.size(); k++) {
    REQUIRE(v[k] == Approx(mvp_control_norms_inf[k]));
  }

  Teuchos::TimeMonitor::summarize();
  return;
}
// ===========================================================================
#if 0
TEST_CASE("hashes for rectangle mesh", "[rectangle]")
{
  testKeo(
      "rectanglesmall",
      4.0,
      1.0,
      {0.25, 2.5, 0.0}
      );
}
#endif
// ============================================================================
TEST_CASE("hashes for pacman mesh", "[pacman]")
{
  testKeo(
      "pacman",
      409.0,
      1.0,
      {4.999111652374270, 5.0, 0.0}
      );
}
// ============================================================================
#if 0
TEST_CASE("hashes for shell mesh", "[shell]")
{
  testKeo(
      "shell",
      5.0,
      1.0,
      {0.5, 0.5, 0.0}
      );
}
#endif
// ============================================================================
TEST_CASE("hashes for sphere mesh", "[sphere]")
{
  testKeo(
      "sphere",
      82.0,
      1.0,
      {0.492403864860535, 0.468303918838501, 0.0}
      );
}
// ============================================================================
#if 0
TEST_CASE("hashes for cube mesh", "[cube]")
{
  testKeo(
      "cubesmall",
      8.0,
      1.0,
      {0.25, 0.25, 0.0}
      );
}
#endif
// ============================================================================
TEST_CASE("hashes for brick mesh", "[brick]")
{
  testKeo(
      "brick-w-hole",
      744.0,
      1.0,
      {2.5, 2.5, 0.0}
      );
}
// ============================================================================
