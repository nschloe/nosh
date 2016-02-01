#ifndef HELPER_HPP
#define HELPER_HPP

#include <boost/any.hpp>
#include <map>

namespace nosh
{
  void
  show_any(const boost::any & any);

  void
  show_map(
    const std::map<std::string, boost::any> & map,
    const int indent = 0
    );
}
#endif
