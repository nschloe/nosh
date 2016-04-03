#ifndef NOSH_MESHREADER_HPP
#define NOSH_MESHREADER_HPP
// =============================================================================
// includes
#include <memory>
#include <set>
#include <string>

#include "mesh.hpp"

namespace nosh
{

std::shared_ptr<nosh::mesh>
read(const std::string & file_name);

} // namespace nosh
// =============================================================================
#endif // NOSH_MESHREADER_HPP
