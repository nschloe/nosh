#include "function.hpp"

void
nosh::
write(
    const nosh::function & x,
    const std::string & filename
    )
{
  //x.mesh->insert_vector(x, "x");
  x.mesh->write(filename);
}
