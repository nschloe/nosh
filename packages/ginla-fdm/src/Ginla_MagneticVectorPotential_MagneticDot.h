#ifndef GINLA_MAGNETICVECTORPOTENTIAL_MAGNETICDOT_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_MAGNETICDOT_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_MagneticVectorPotential_Virtual.h"

namespace Ginla {
  namespace MagneticVectorPotential {

class MagneticDot:
  public Virtual
{
public:
  MagneticDot( double mu );

  virtual
  ~MagneticDot();

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

  virtual
  double
  getAx(const Point & x) const;

  virtual
  double
  getAy(const Point & x) const;

  virtual
  double
  getAz(const Point & x) const;

  virtual
  Teuchos::RCP<Point>
  getDADMu(const Point & x ) const;

  virtual
  double
  getDAxDMu(const Point & x ) const;

  virtual
  double
  getDAyDMu(const Point & x ) const;

  virtual
  double
  getDAzDMu(const Point & x ) const;

protected:
private:
  //! Radius of the magnetic dot.
  const double magnetRadius_;

  //! \f$z\f$-distance between the SC sample and lower surface of magnetic disk.
  const double zz1_;

  //! \f$z\f$-distance between the SC sample and upper surface of magnetic disk.
  const double zz2_;

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_MAGNETICDOT_H_
