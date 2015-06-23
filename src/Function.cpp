
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
  try {
    x.mesh->write();
  }
  catch (const std::runtime_error& error)
  {
    // Ignore std::runtime_error.
    // Those are hopefully restricted to
    // ```
    //[ex_put_time] Error: failed to store time value in file id 65536
    //    exerrval = -130
    //terminate called after throwing an instance of 'std::runtime_error'
    //  what():  Exodus error (-130)NetCDF: Attempt to extend dataset during NC_INDEPENDENT I/O operation. Use nc_var_par_access to set mode NC_COLLECTIVE before extending variable. at line 853 in file 'Ioex_DatabaseIO.C 2015/04/13' Please report to gdsjaar@sandia.gov if you need help.
    // ```
    // Greg should be able to handle those soon.
  }
  auto mesh2 = Nosh::read(tmpFile, {"x"});
  mesh2->insertVector(x, "x");
  mesh2->openFile(fileName);
  try {
    mesh2->write();
  }
  catch (const std::runtime_error& error)
  {
    // see above
  }
}
