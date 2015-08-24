// @HEADER
//
//    Query routines for a vector potential with explicitly given values.
//    Copyright (C) 2011, 2012  Nico Schlömer
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
#ifndef NOSH_VECTORFIELD_EXPLICITVALUES_H_
#define NOSH_VECTORFIELD_EXPLICITVALUES_H_

#include <map>
#include <string>

#include "vector_field_base.hpp"
#include "mesh.hpp"

namespace nosh
{
namespace vector_field
{
class explicit_values : public base
{
public:
  explicit_values(
      const nosh::mesh &mesh,
      const std::string &field_name,
      const double mu
      );

  virtual
  ~explicit_values();

  virtual
  void
  set_parameters(const std::map<std::string, double> & params);

  //! Get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  get_parameters() const;

  double
  get_edge_projection(const unsigned int edge_index) const;

  double
  get_d_edge_projection_dp(
      const unsigned int edge_index,
      const std::string & dParamName
      ) const;

  virtual
  Eigen::Vector3d
  eval(const Eigen::Vector3d & x) const {
    // mu * 0.5  * X.cross(B)
    return mu_ * 0.5 * Eigen::Vector3d({-x[1], x[0], 0.0});
  };

  virtual
  unsigned int
  degree() const {
    return 1;
  };

protected:
private:
  double mu_;

  std::vector<double> edgeProjectionCache_;
};
} // namespace vector_field
} // namespace nosh
#endif // NOSH_VECTORFIELD_EXPLICITVALUES_H_
