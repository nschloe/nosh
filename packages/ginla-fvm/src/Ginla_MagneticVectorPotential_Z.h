#ifndef GINLA_MAGNETICVECTORPOTENTIAL_Z_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_Z_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_MagneticVectorPotential_Virtual.h"

namespace Ginla {
  namespace MagneticVectorPotential {

class Z:
  public Virtual
{
public:
  Z( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
     double mu );

  virtual
  ~Z();

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
  getRawA_(const Point & x) const;

  void
  initializeEdgeMidpointProjectionCache_() const;

private:

  double mu_;

  const Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeMidpointProjectionCache_;
  mutable bool edgeMidpointProjectionCacheUpToDate_;

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_Z_H_
