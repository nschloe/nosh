#ifndef NOSH_CSVWRITER_H
#define NOSH_CSVWRITER_H

#include <string>

#include <Teuchos_ParameterList.hpp>

namespace nosh
{

class csv_writer
{
public:
  //! Default constructor.
  csv_writer(
      const std::string &file_name,
      const std::string &delimeter = ","
      );

  //! Destructor.
  virtual
  ~csv_writer();

  void
  write_header(const Teuchos::ParameterList & pList) const;

  void
  write_row(const Teuchos::ParameterList & pList) const;

protected:
private:
  //! File stream for the statistics.
  mutable std::ofstream fileStream_;

  const std::string delimeter_;
  const std::string headerStart_;

  const unsigned int doublePrec_;
  const unsigned int doubleColumnWidth_;
  const unsigned int intColumnWidth_;
};
} // namespace nosh

#endif // NOSH_CSVWRITER_H
