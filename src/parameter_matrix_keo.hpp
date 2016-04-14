#ifndef NOSH_MATRIXBUILDER_KEO_H
#define NOSH_MATRIXBUILDER_KEO_H

#include <map>
#include <string>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include <Tpetra_CrsMatrix.hpp>

#include <Eigen/Dense>

#include "mesh.hpp"
#include "parameter_object.hpp"

// forward declarations
namespace nosh
{
class mesh;
namespace scalar_field
{
class base;
}
namespace vector_field
{
class base;
}
} // namespace nosh

namespace nosh
{
namespace parameter_matrix
{

class keo: public nosh::parameter_object, public Tpetra::CrsMatrix<double,int,int>
{
public:
  keo(
      const std::shared_ptr<const nosh::mesh> &mesh,
      const std::shared_ptr<const nosh::scalar_field::base> &thickness,
      const std::shared_ptr<nosh::vector_field::base> &mvp
     );

  // Destructor.
  ~keo();

  //! Gets the initial parameters from this module.
  virtual
  const std::map<std::string, double>
  get_scalar_parameters() const;

protected:
private:
  void
  refill_(
      const std::map<std::string, double> & scalar_params,
      const std::map<std::string, std::shared_ptr<const Tpetra::Vector<double, int, int>>> & vector_params
      );

  void
  build_alpha_cache_(
      const std::vector<edge> & edges,
      const std::vector<nosh::mesh::edge_data> & edge_data
      ) const;

  double
  integrate1d_(
      const nosh::vector_field::base & f,
      const Eigen::Vector3d & x0,
      const Eigen::Vector3d & x1
      ) const;

private:
  const std::shared_ptr<const nosh::mesh> mesh_;
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> keo_fill_time_;
#endif
  const std::shared_ptr<const nosh::scalar_field::base> thickness_;
  const std::shared_ptr<nosh::vector_field::base> mvp_;

  mutable std::vector<double> alpha_cache_;
  mutable bool alpha_cache_up_to_date_;
};
} // namespace parameter_matrix
} // namespace nosh

#endif // NOSH_MATRIXBUILDER_KEO_H
