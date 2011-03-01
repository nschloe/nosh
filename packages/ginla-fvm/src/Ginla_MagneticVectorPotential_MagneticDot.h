#ifndef GINLA_MAGNETICVECTORPOTENTIAL_MAGNETICDOT_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_MAGNETICDOT_H_

#include "Ginla_MagneticVectorPotential_Virtual.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

namespace Ginla {
  namespace MagneticVectorPotential {

class MagneticDot:
  public Virtual
{
public:
  MagneticDot( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
               double magnetRadius,
               double mu
             );

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
  getAEdgeMidpointProjection( const unsigned int cellIndex,
                              const unsigned int edgeIndex
                            ) const;

protected:

private:

  Teuchos::RCP<Point>
  getRawA_( const Point & x ) const;

  void
  initializeEdgeMidpointProjectionCache_() const;

private:

  double mu_;

  //! Radius of the magnetic dot.
  const double magnetRadius_;

  //! \f$z\f$-distance between the SC sample and lower surface of magnetic disk.
  const double zz1_;

  //! \f$z\f$-distance between the SC sample and upper surface of magnetic disk.
  const double zz2_;

  const Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeMidpointProjectionCache_;
  mutable bool edgeMidpointProjectionCacheUpToDate_;

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_MAGNETICDOT_H_
