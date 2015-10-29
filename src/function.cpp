
#include "function.hpp"

#include "mesh_reader.hpp"

void
nosh::
write(
    const nosh::function & x,
    const std::string & file_name
    )
{
  // 2015-06-20:
  // STK has to know the fields that are written out with a mesh *before* the
  // mesh is read. Hopefully, the reader was smart enough allocate some
  // space...
  const std::set<std::string> names = x.mesh->allocated_vector_names;
  if (!names.empty()) {
    std::string allocated_vector_name;
    for (const auto& name: names) {
      // any name will do
      allocated_vector_name = name;
      break;
    }
    x.mesh->insert_vector(x, allocated_vector_name);
    x.mesh->open_file(file_name);
    x.mesh->write();
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "STK meshes need to have space reserved for writing at read time."
        );
  }
}
