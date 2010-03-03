/*
 * MagneticVectorPotential.h
 *
 *  Created on: Nov 5, 2009
 *      Author: nico
 */

#ifndef MAGNETICVECTORPOTENTIAL_H_
#define MAGNETICVECTORPOTENTIAL_H_

#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

namespace GL {
  namespace MagneticVectorPotential {

class Centered
{
public:
  // TODO Remove edge length from this class.
  Centered( double h0,
            double edgeLength );

  virtual
  ~Centered();

  void
  setH0( const double h0 );

  void
  setEdgeLength( const double edgeLength );

  double
  getH0() const;

  void
  setParameters( const LOCA::ParameterVector & p );

  Teuchos::RCP<Teuchos::Array<double> >
  getA(const Teuchos::Array<double> & x ) const;

  double
  getAx(const Teuchos::Array<double> & x) const;

  double
  getAy(const Teuchos::Array<double> & x) const;

protected:

private:

  double h0_;
  double edgeLength_; // TODO: Remove geometric info from magnetic vector potential

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif /* MAGNETICVECTORPOTENTIAL_H_ */
