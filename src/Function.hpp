
#include <memory>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "Mesh.hpp"

namespace Nosh {

  class Function:
    public Tpetra::Vector<double,int,int>
  {
    public:
      Function(const std::shared_ptr<Nosh::Mesh> inMesh):
        Tpetra::Vector<double,int,int>(Teuchos::rcp(inMesh->getMap())),
        mesh(inMesh)
      {
      };

      virtual
      ~Function() {};

    public:
      const std::shared_ptr<Nosh::Mesh> mesh;
  };

  // Helper functions
  void
  write(
      const Nosh::Function & x,
      const std::string & fileName
      );

}
