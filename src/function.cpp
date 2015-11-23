
#include "function.hpp"

#include "mesh_reader.hpp"

void
nosh::
write(
    const nosh::function & x,
    const std::string & filename
    )
{
  x.mesh->write(filename);
}
