#ifndef GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <LOCA_Parameter_Vector.H>

// typedef just like in ...
typedef Teuchos::Tuple<double,3> Point;

namespace Ginla {
  namespace MagneticVectorPotential {

class Virtual
{
public:
  Virtual( double mu );

  virtual
  ~Virtual();

  bool
  setMu( const double mu );
  
  //! Sets the parameters in this module.
  //! @return Indicates whether the internal values have changed.
  virtual
  bool
  setParameters( const LOCA::ParameterVector & p ) = 0;

  virtual
  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const = 0;

  virtual
  Teuchos::RCP<Point>
  getA(const Point & x ) const = 0;

  virtual
  double
  getAx(const Point & x) const = 0;

  virtual
  double
  getAy(const Point & x) const = 0;
  
  virtual
  double
  getAz(const Point & x) const = 0;
  
  virtual
  Teuchos::RCP<Point>
  getDADMu(const Point & x ) const = 0;

  virtual
  double
  getDAxDMu(const Point & x ) const = 0;

  virtual
  double
  getDAyDMu(const Point & x ) const = 0;
  
  virtual
  double
  getDAzDMu(const Point & x ) const = 0;

protected:
  double mu_;

private:
};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
