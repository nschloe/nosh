#ifndef HELPER_HPP
#define HELPER_HPP

#include <iostream>
#include <map>
#include <string>

namespace nosh
{
  void
  show_map(
    const std::map<std::string, boost::any> & map,
    const int indent = 0
    )
  {
    std::cout << "{\n";
    for (const auto &p: map) {
      std::cout << std::string(indent + 2, ' ') << p.first << ": ";
      if (p.second.type() == typeid(int)) {
        std::cout << boost::any_cast<int>(p.second);
      } else if (p.second.type() == typeid(double)) {
        std::cout << boost::any_cast<double>(p.second);
      } else if (p.second.type() == typeid(std::string)) {
        std::cout << boost::any_cast<std::string>(p.second);
      } else if (p.second.type() == typeid(const char*)) {
        std::cout << boost::any_cast<const char*>(p.second);
      } else if (p.second.type() == typeid(std::map<std::string, boost::any>)) {
        show_map(
            boost::any_cast<std::map<std::string, boost::any>>(p.second),
            indent + 2
            );
      } else {
        std::cout << "[unhandled type]";
      }
      std::cout << std::endl;
    }
    std::cout << std::string(indent, ' ') << "}";
  }
}
#endif
