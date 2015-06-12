
#include <Tpetra_Vector.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "Mesh.hpp"

namespace Nosh {

  class Function:
    public Tpetra::Vector<double,int,int>
  {
    public:
      Function(const std::shared_ptr<const Nosh::Mesh> & mesh):
        Tpetra::Vector<double,int,int>(Teuchos::rcp(mesh->getMap())),
        mesh_(mesh)
      {
      };

      virtual
      ~Function() {};

    private:
      const std::shared_ptr<const Nosh::Mesh> mesh_;
  };

}
