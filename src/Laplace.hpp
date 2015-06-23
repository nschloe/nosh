// @HEADER
//
//    Builder class for the Laplace operator.
//    Copyright (C) 2012  Nico Schlömer
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
#ifndef NOSH_LAPLACE_H
#define NOSH_LAPLACE_H

#include "EdgeOperator.hpp"

namespace Nosh
{
  class Laplace:
    public EdgeOperator
  {
    public:
      Laplace(
          const std::shared_ptr<const Nosh::Mesh> & _mesh,
          const std::set<std::shared_ptr<const Nosh::DirichletBC>> & _bcs
          ):
        EdgeOperator(_mesh, _bcs)
      {
        this->fill_();
      };

      virtual
      ~Laplace()
      {};

    protected:
      virtual
        std::vector<std::vector<double>>
        edgeContrib(
            const double edgeCoefficient,
            const double controlVolume0,
            const double controlVolume1
            ) const
        {
          (void) controlVolume0;
          (void) controlVolume1;
          return {
            { edgeCoefficient, -edgeCoefficient},
            {-edgeCoefficient,  edgeCoefficient}
          };
        }

  };
} // namespace Nosh

#endif // NOSH_LAPLACE_H
