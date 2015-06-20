
#include "Function.hpp"

#include "MeshReader.hpp"

void
Nosh::
write(
    const Nosh::Function & x,
    const std::string & fileName
    )
{
  // 2015-06-20:
  // STK has to know the fields that are written out with a mesh *before* the
  // mesh is read. This makes writing any function out with a mesh impossible.
  // Well, the workaround is to write out the mesh to a temporary file, read it
  // in, claiming there is a field "x", and then writing it out again.
  const std::string tmpFile = "/tmp/mesh.e";
  x.mesh->openFile(tmpFile);
  x.mesh->write();
  auto mesh2 = Nosh::read(tmpFile, {"x"});
  mesh2->insertVector(x, "x");
  mesh2->openFile(fileName);
  mesh2->write();
}
