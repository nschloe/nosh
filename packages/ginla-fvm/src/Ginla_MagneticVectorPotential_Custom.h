#ifndef GINLA_MAGNETICVECTORPOTENTIAL_CUSTOM_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_CUSTOM_H_

#include "Ginla_MagneticVectorPotential_Virtual.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Epetra_MultiVector.h>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

namespace Ginla {
  namespace MagneticVectorPotential {

class Custom:
  public Virtual
{
public:
  Custom( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
          const Teuchos::RCP<const Epetra_MultiVector>  & mvp,
          double mu
        );

  virtual
  ~Custom();

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

  void
  initializeEdgeMidpointProjectionCache_() const;

private:

  double mu_;

  const Teuchos::RCP<const Epetra_MultiVector> mvp_;

  const Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeMidpointProjectionCache_;
  mutable bool edgeMidpointProjectionCacheUpToDate_;

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_CUSTOM_H_
