#include <catch.hpp>

#include <string>

#include <Teuchos_DefaultComm.hpp>

#include <nosh.hpp>

// =============================================================================
void
testMesh(
    const std::string & input_filename_base,
    const unsigned int control_num_nodes,
    const double control_vol_norm_1,
    const double control_vol_norm_2,
    const double control_vol_norm_inf
    )
{
  auto comm =  Teuchos::DefaultComm<int>::getComm();
  const int size = comm->getSize();
  const std::string input_filename = (size == 1) ?
    "data/" + input_filename_base + ".h5m" :
    "data/" + input_filename_base + "-" + std::to_string(size) + ".h5m"
    ;

  auto mesh = nosh::read(input_filename);

  const unsigned int num_nodes = mesh->map()->getGlobalNumElements();
  REQUIRE(num_nodes == control_num_nodes);

  const auto control_vols = mesh->control_volumes();
  REQUIRE(control_vols->norm1() == Approx(control_vol_norm_1));
  REQUIRE(control_vols->norm2() == Approx(control_vol_norm_2));
  REQUIRE(control_vols->normInf() == Approx(control_vol_norm_inf));

  return;
}
// ===========================================================================
#if 0
TEST_CASE("hashes for rectangle mesh", "[rectangle]")
{
  testMesh(
      "rectanglesmall",
      4,
      10.0,
      5.0,
      2.5
      );
}
#endif
// ============================================================================
TEST_CASE("hashes for pacman mesh", "[pacman]")
{
  testMesh(
      "pacman",
      409,
      302.5227007210103,
      15.38575790933914,
      1.127797467043659
      );
}
// ============================================================================
#if 0
TEST_CASE("hashes for shell mesh", "[shell]")
{
  testMesh(
      "shell",
      5,
      3.46410161513775,
      1.63299316185545,
      1.15470053837925
      );
}
#endif
// ============================================================================
TEST_CASE("hashes for sphere mesh", "[sphere]")
{
  testMesh(
      "sphere",
      82,
      11.9741927059035,
      1.39047542328083,
      0.198927169088121
      );
}
// ===========================================================================
#if 0
TEST_CASE("hashes for cube mesh", "[cube]")
{
  testMesh(
      "cubesmall",
      8,
      10.0,
      3.535533905932738,
      1.25
      );
}
#endif
// ============================================================================
TEST_CASE("hashes for brick mesh", "[brick]")
{
  testMesh(
      "brick-w-hole",
      744,
      388.686291694641,
      16.6614019419857,
      1.46847345474977
      );
}
// ============================================================================
