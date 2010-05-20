/*
 * MagneticVectorPotential.h
 *
 *  Created on: Nov 5, 2009
 *      Author: nico
 */

#ifndef GINLA_MAGNETICVECTORPOTENTIAL_CENTERED_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_CENTERED_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <LOCA_Parameter_Vector.H>

// typedef just like in Recti
typedef Teuchos::Tuple<double,2> DoubleTuple;


namespace Ginla {
  namespace MagneticVectorPotential {

class Centered
{
public:
  // TODO Remove edge length from this class.
  Centered( double h0,
            double edgeLength );

  virtual
  ~Centered();

  //! Sets the parameters in this module.
  //! @return Indicates whether the internal values have changed.
  bool
  setParameters( const LOCA::ParameterVector & p );
  
  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const;

  Teuchos::RCP<DoubleTuple>
  getA(const DoubleTuple & x ) const;

  double
  getAx(const DoubleTuple & x) const;

  double
  getAy(const DoubleTuple & x) const;
  
  Teuchos::RCP<DoubleTuple>
  getDADh0(const DoubleTuple & x ) const;

  double
  getDAxDh0(const DoubleTuple & x ) const;

  double
  getDAyDh0(const DoubleTuple & x ) const;

protected:

private:

  double h0_;
  double edgeLength_; // TODO: Remove geometric info from magnetic vector potential

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif /* GINLA_MAGNETICVECTORPOTENTIAL_CENTERED_H_ */
