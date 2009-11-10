/*
 * MagneticVectorPotential.h
 *
 *  Created on: Nov 5, 2009
 *      Author: nico
 */

#ifndef MAGNETICVECTORPOTENTIAL_H_
#define MAGNETICVECTORPOTENTIAL_H_

#include <Teuchos_Array.hpp>

class MagneticVectorPotential
{
public:
  // TODO:
  // Remove edge length from this class.
  MagneticVectorPotential( double h0,
                           double edgeLength );

  virtual
  ~MagneticVectorPotential();

  void
  setH0( const double h0 );

  double
  getH0() const;

  Teuchos::RCP<Teuchos::Array<double> >
  getA(const Teuchos::Array<double> & x ) const;

  double
  getAx(const Teuchos::Array<double> & x) const;

  double
  getAy(const Teuchos::Array<double> & x) const;

protected:

private:

  double h0_;
  double edgeLength_; // TODO: remove this

};

#endif /* MAGNETICVECTORPOTENTIAL_H_ */
