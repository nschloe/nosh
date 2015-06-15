
#include "Function.hpp"

void
Nosh::
write(
    const Nosh::Function & x,
    const std::string & fileName
    )
{
  x.mesh->insertVector(x, "x");
  x.mesh->openFile(fileName);
  x.mesh->write();
}
