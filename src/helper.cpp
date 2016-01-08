#include "helper.hpp"

#include <iostream>
#include <map>
#include <string>

namespace nosh
{
  void
  show_any(const boost::any & any)
  {
    if (any.type() == typeid(int)) {
      std::cout << boost::any_cast<int>(any);
    } else if (any.type() == typeid(double)) {
      std::cout << boost::any_cast<double>(any);
    } else if (any.type() == typeid(std::string)) {
      std::cout << boost::any_cast<std::string>(any);
    } else if (any.type() == typeid(const char*)) {
      std::cout << boost::any_cast<const char*>(any);
    } else if (any.type() == typeid(std::map<std::string, boost::any>)) {
      show_map(
          boost::any_cast<std::map<std::string, boost::any>>(any),
          2
          );
    } else {
      std::cout << "[unhandled type]";
    }
  }

  void
  show_map(
    const std::map<std::string, boost::any> & map,
    const int indent
    )
  {
    std::cout << "{\n";
    for (const auto &p: map) {
      std::cout << std::string(indent + 2, ' ') << p.first << ": ";
      show_any(p.second);
      std::cout << std::endl;
    }
    std::cout << std::string(indent, ' ') << "}";
  }
}
