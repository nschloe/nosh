#ifndef GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <LOCA_Parameter_Vector.H>

// forward declarations
namespace Ginla {
  namespace EpetraFVM {
    class StkMesh;
  }
}

typedef Teuchos::Tuple<double,3> Point;

namespace Ginla {
  namespace MagneticVectorPotential {

class Virtual
{
public:
  Virtual( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh );

  virtual
  ~Virtual();

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

  //! Return the projection of the magnetic vector potential onto the
  //! edge with index \c edgeIndex at the midpoint of the edge.
  virtual
  double
  getAEdgeMidpointProjection( const unsigned int cellIndex,
                              const unsigned int edgeIndex
                            ) const = 0;

protected:
  const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> mesh_;
private:
};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
