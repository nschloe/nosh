#include "jacobian_operator.hpp"

#include <map>
#include <string>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "parameter_matrix_keo.hpp"
#include "mesh.hpp"
#include "scalar_field_base.hpp"

namespace nosh
{
// =============================================================================
jacobian_operator::
jacobian_operator(
    const std::shared_ptr<const nosh::mesh> &mesh,
    const std::shared_ptr<const nosh::scalar_field::base> &scalar_potential,
    const std::shared_ptr<const nosh::scalar_field::base> &thickness,
    const std::shared_ptr<nosh::parameter_matrix::keo> &keo
    ) :
  mesh_(mesh),
  scalar_potential_(scalar_potential),
  thickness_(thickness),
  keo_(keo),
  diag0_(Teuchos::rcp(mesh->complex_map())),
  diag1b_(mesh->control_volumes()->getMap())
{
}
// =============================================================================
jacobian_operator::
~jacobian_operator()
{
}
// =============================================================================
void
jacobian_operator::
apply(
    const Tpetra::MultiVector<double,int,int> &X,
    Tpetra::MultiVector<double,int,int> &Y,
    Teuchos::ETransp mode,
    double alpha,
    double beta
    ) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      mode != Teuchos::NO_TRANS,
      "Only untransposed applies supported."
      );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      alpha != 1.0,
      "Only alpha==1.0 supported."
      );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      beta != 0.0,
      "Only beta==0.0 supported."
      )
  // Add the terms corresponding to the nonlinear terms.
  // A = K + I * thickness * (V + g * 2*|psi|^2)
  // B = g * diag(thickness * psi^2)

  // Y = K*X
  keo_->apply(X, Y);

  const int num_my_points = mesh_->control_volumes()->getLocalLength();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*num_my_points, X.getLocalLength());
#endif

  auto d0_data = diag0_.getData();
  auto d1b_data = diag1b_.getData();

#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*num_my_points, d0_data.size());
  TEUCHOS_ASSERT_EQUALITY(num_my_points, d1b_data.size());
#endif

  for (std::size_t i = 0; i < X.getNumVectors(); i++) {
    auto x_data = X.getVector(i)->getData();
    auto y_data = Y.getVectorNonConst(i)->getDataNonConst();
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(2*num_my_points, x_data.size());
    TEUCHOS_ASSERT_EQUALITY(2*num_my_points, y_data.size());
#endif
    // For the parts Re(psi)Im(phi), Im(psi)Re(phi), the (2*k+1)th component of
    // X needs to be summed into the (2k)th component of Y, likewise for (2k)
    // -> (2k+1).
    // The Epetra class cannot currently handle this situation (e.g., by
    // Tpetra::Vector<double,int,int>::Multiply()), so we need to access the
    // vector entries one-by-one. And then, while we're at it, let's include
    // all the other terms in the loop too. (It would actually be possible to
    // have the terms 2k/2k and 2k+1/2k+1 handled by Multiply().
    for (int k = 0; k < num_my_points; k++) {
      y_data[2*k] += d0_data[2*k] * x_data[2*k]
                  + d1b_data[k]  * x_data[2*k+1];
      y_data[2*k+1] += d1b_data[k]    * x_data[2*k]
                    + d0_data[2*k+1] * x_data[2*k+1];
    }
  }

//    // take care of the shifting
//    if (alpha_ != 0.0 || beta_ != -1.0)
//    TEUCHOS_ASSERT_EQUALITY(0, Y.Update(alpha_, X, -beta_));

  return;
}
// =============================================================================
void
jacobian_operator::
rebuild(
    const std::map<std::string, double> & params,
    const Tpetra::Vector<double,int,int> & current_x
    )
{
  // Fill the KEO.
  // It is certainly a debatable design decision to have our own KEO in
  // jacobian_operator and not live of the cache of the builder. On the one
  // hand, in a typical continuation context, the same matrix is used in
  // compute_f, the preconditioner, and here. It would then be sufficient to
  // store the matrix at one common place (e.g., the builder).  This might
  // however lead to complications in a situation like the following:
  //
  //   1. The matrix is requested by the Jacobian operator, a pointer to the
  //      common storage place is returned.
  //   2. The matrix is requested by compute_f() with different parameters.
  //      Now also the Jacobian operator's instance has the altered matrix.
  //   3. The Jacobian operator uses the matrix.
  //
  // One might argue that this situation is does not occur in the given
  // context, but really the code shouldn't make any assumptions about it.
  // Besides, the matrix copy that happens in fill is not of much concern
  // computationally. Should this ever become an issue, revisit.
  keo_->set_parameters(params, {});

  // Rebuild diagonals.
  this->rebuild_diags_(params, current_x);

  return;
}
// =============================================================================
void
jacobian_operator::
rebuild_diags_(
    const std::map<std::string, double> & params,
    const Tpetra::Vector<double,int,int>  &x
    )
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(scalar_potential_);
#endif

  const auto & control_volumes = *(mesh_->control_volumes());

  const double g = params.at("g");

  const auto thickness_values = thickness_->get_v(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(control_volumes.getMap()->isSameAs(*thickness_values.getMap()));
#endif

  const auto scalar_potential_values = scalar_potential_->get_v(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(
      control_volumes.getMap()->isSameAs(*scalar_potential_values.getMap())
      );
#endif

  auto x_data = x.getData();
  auto c_data = control_volumes.getData();
  auto t_data = thickness_values.getData();
  auto s_data = scalar_potential_values.getData();

  auto d0_data = diag0_.getDataNonConst();
  auto d1b_data = diag1b_.getDataNonConst();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(c_data.size(), t_data.size());
  TEUCHOS_ASSERT_EQUALITY(t_data.size(), s_data.size());
  TEUCHOS_ASSERT_EQUALITY(2*s_data.size(), x_data.size());
  TEUCHOS_ASSERT_EQUALITY(s_data.size(), d1b_data.size());
#endif

  for (decltype(c_data)::size_type k = 0; k < c_data.size(); k++) {
    // rebuild diag0
    const double alpha = c_data[k] * t_data[k]
      * (s_data[k]
          + g * 2.0 * (x_data[2*k]*x_data[2*k] + x_data[2*k+1]*x_data[2*k+1])
        );
    const double realX2 = g * c_data[k] * t_data[k]
      * (x_data[2*k]*x_data[2*k] - x_data[2*k+1]*x_data[2*k+1]);
    d0_data[2*k]   = alpha + realX2;
    d0_data[2*k+1] = alpha - realX2;

    // rebuild diag1b
    d1b_data[k] = g * c_data[k] * t_data[k] * (2.0 * x_data[2*k] * x_data[2*k+1]);
  }

  return;
}
// =============================================================================
} // namespace nosh
