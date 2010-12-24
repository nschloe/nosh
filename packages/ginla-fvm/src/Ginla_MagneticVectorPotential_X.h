#ifndef GINLA_MAGNETICVECTORPOTENTIAL_X_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_X_H_

#include "Ginla_MagneticVectorPotential_Virtual.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <LOCA_Parameter_Vector.H>

namespace Ginla {
  namespace MagneticVectorPotential {

class X:
  public Virtual
{
public:
  X( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
     double mu );

  virtual
  ~X();

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

  double mu_;

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_X_H_
