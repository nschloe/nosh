#ifndef NOSH_FUNCTION_HPP
#define NOSH_FUNCTION_HPP

#include <memory>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "mesh.hpp"

namespace nosh {

  class function:
    public Tpetra::Vector<double,int,int>
  {
    public:
      function(const std::shared_ptr<nosh::mesh> in_mesh):
        Tpetra::Vector<double,int,int>(Teuchos::rcp(in_mesh->map())),
        mesh(in_mesh)
      {
      };

      virtual
      ~function() {};

    public:
      const std::shared_ptr<nosh::mesh> mesh;
  };

  // Helper functions
  void
  write(
      const nosh::function & x,
      const std::string & file_name
      );

}
#endif // NOSH_FUNCTION_HPP
