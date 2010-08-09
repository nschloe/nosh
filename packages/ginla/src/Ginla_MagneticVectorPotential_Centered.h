#ifndef GINLA_MAGNETICVECTORPOTENTIAL_CENTERED_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_CENTERED_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <LOCA_Parameter_Vector.H>

// typedef just like in ...
typedef Teuchos::Tuple<double,3> Point;

namespace Ginla {
  namespace MagneticVectorPotential {

class Centered
{
public:
  Centered( double h0 );

  virtual
  ~Centered();

  //! Sets the parameters in this module.
  //! @return Indicates whether the internal values have changed.
  bool
  setParameters( const LOCA::ParameterVector & p );
  
  bool
  setH0( const double mu );
  
  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const;

  Teuchos::RCP<Point>
  getA(const Point & x ) const;

  double
  getAx(const Point & x) const;

  double
  getAy(const Point & x) const;
  
  double
  getAz(const Point & x) const;
  
  Teuchos::RCP<Point>
  getDADh0(const Point & x ) const;

  double
  getDAxDh0(const Point & x ) const;

  double
  getDAyDh0(const Point & x ) const;
  
  double
  getDAzDh0(const Point & x ) const;

protected:

private:

  double h0_;

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif /* GINLA_MAGNETICVECTORPOTENTIAL_CENTERED_H_ */
