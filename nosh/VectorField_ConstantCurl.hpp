// @HEADER
//
//    Query routines for the vector potential associated with a constant curl field.
//    Copyright (C) 2011, 2012  Nico Schl\"omer
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
#ifndef NOSH_VECTORFIELD_CONSTANTCURL_H_
#define NOSH_VECTORFIELD_CONSTANTCURL_H_

#include <map>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Epetra_MultiVector.h>
#include <Teuchos_Array.hpp>

#include "nosh/VectorField_Virtual.hpp"
#include "nosh/StkMesh.hpp"

typedef Teuchos::SerialDenseVector<int, double> DoubleVector;

namespace Nosh
{
namespace VectorField
{
class ConstantCurl : public Virtual
{
public:
  ConstantCurl(const Teuchos::RCP<Nosh::StkMesh> &mesh,
               const Teuchos::RCP<DoubleVector> &b,
               const Teuchos::RCP<DoubleVector> &u = Teuchos::null
              );

  virtual
  ~ConstantCurl();

//! Get the parameter names and initial values.
  virtual
  const std::map<std::string, double>
  getInitialParameters() const;

  virtual
  double
  getEdgeProjection(const unsigned int edgeIndex,
                    const std::map<std::string, double> & params
                  ) const;

  virtual
  double
  getDEdgeProjectionDp(const unsigned int edgeIndex,
                       const std::map<std::string, double> & params,
                       const std::string & dParam
                     ) const;

protected:
private:
  DoubleVector
  getRawA_(const DoubleVector &x) const;

  DoubleVector
  getRawDADTheta_(const DoubleVector &x) const;

  void
  initializeEdgeCache_() const;

  void
  rotate_(DoubleVector &v,
          const DoubleVector &u,
          const double theta
        ) const;

  void
  dRotateDTheta_(DoubleVector &v,
                 const DoubleVector &u,
                 const double theta
               ) const;

  DoubleVector
  crossProduct_(const DoubleVector u,
                const DoubleVector v
              ) const;

private:
  const Teuchos::RCP<Nosh::StkMesh> mesh_;
  const Teuchos::RCP<const DoubleVector> b_;
  const Teuchos::RCP<const DoubleVector> u_;
  mutable DoubleVector rotatedBCache_;
  mutable double rotatedBCacheAngle_;
  mutable DoubleVector dRotatedBDThetaCache_;
  mutable double rotateddBdThetaCacheAngle_;

  Teuchos::ArrayRCP<DoubleVector> edgeCache_;
  mutable bool edgeCacheUptodate_;
};
} // namespace VectorField
} // namespace Nosh
#endif // NOSH_VECTORFIELD_CONSTANTCURL_H_
