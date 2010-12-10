#ifndef GINLA_MAGNETICVECTORPOTENTIAL_X_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_X_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_MagneticVectorPotential_Virtual.h"

namespace Ginla {
  namespace MagneticVectorPotential {

class Spherical:
  public Virtual
{
public:
  Spherical( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
             double mu,
             double phi,
             double theta );

  virtual
  ~Spherical();

  //! Sets the parameters in this module.
  //! @return Indicates whether the internal values have changed.
  virtual
  bool
  setParameters( const LOCA::ParameterVector & p );

  virtual
  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const;

  virtual
  Teuchos::RCP<Point>
  getA(const Point & x ) const;

  double
  getAEdgeMidpointProjection( const unsigned int cellIndex,
                              const unsigned int edgeIndex
                            ) const;

protected:
private:

  Teuchos::RCP<Point>
  getRawA_(const Point & x) const;

  void
  initializeEdgeMidpointProjectionCache_() const;

  double
  getAx(const Point & x) const;

  double
  getAy(const Point & x) const;

  double
  getAz(const Point & x) const;

  Teuchos::RCP<Point>
  getDADMu(const Point & x) const;

  double
  getDAxDMu(const Point & x) const;

  double
  getDAyDMu(const Point & x) const;

  double
  getDAzDMu(const Point & x) const;

private:

  double mu_;
  double phi_;
  double theta_;

  const Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeMidpointProjectionCache_;
  mutable bool edgeMidpointProjectionCacheUpToDate_;

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_X_H_
