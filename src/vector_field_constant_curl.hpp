#ifndef NOSH_VECTORFIELD_CONSTANTCURL_H_
#define NOSH_VECTORFIELD_CONSTANTCURL_H_

#include <map>
#include <string>

#include <Teuchos_RCP.hpp>

#include "vector_field_base.hpp"
#include "mesh.hpp"

namespace nosh
{
namespace vector_field
{
class constantCurl : public base
{
public:
  constantCurl(std::shared_ptr<nosh::mesh> mesh,
               const std::shared_ptr<Eigen::Vector3d> &b,
               const std::shared_ptr<Eigen::Vector3d> &u = nullptr
              );

  ~constantCurl() override;

  void
  set_parameters(const std::map<std::string, double> & params) override;

  //! get parameter names and initial values.
  const std::map<std::string, double>
  get_scalar_parameters() const override;

  double
  get_edge_projection(const unsigned int edge_index) const override;

  double
  get_d_edge_projection_dp(
      const unsigned int edge_index,
      const std::string & param_name
      ) const override;

protected:
private:
  Eigen::Vector3d
  getRawA_(const Eigen::Vector3d &x) const;

  Eigen::Vector3d
  getRawDADTheta_(const Eigen::Vector3d &x) const;

  void
  initializeEdgeCache_() const;

  void
    rotate_(
        Eigen::Vector3d &v,
        const Eigen::Vector3d &u,
        const double theta
        ) const;

  void
    dRotateDTheta_(
        Eigen::Vector3d &v,
        const Eigen::Vector3d &u,
        const double theta
        ) const;

private:
  const std::shared_ptr<nosh::mesh> mesh_;
  const std::shared_ptr<const Eigen::Vector3d> b_;
  const std::shared_ptr<const Eigen::Vector3d> u_;
  mutable Eigen::Vector3d rotatedBCache_;
  mutable double rotatedBCacheAngle_;
  mutable Eigen::Vector3d dRotatedBDThetaCache_;
  mutable double rotateddBdThetaCacheAngle_;

  mutable std::vector<Eigen::Vector3d> edgeCache_;
  mutable bool edgeCacheUptodate_;

  double mu_;
  double theta_;
};
} // namespace vector_field
} // namespace nosh
#endif // NOSH_VECTORFIELD_CONSTANTCURL_H_
