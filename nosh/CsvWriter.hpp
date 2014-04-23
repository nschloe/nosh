// @HEADER
//
//    Helper class for writing out statistics.
//    Copyright (C) 2010--2012  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
#ifndef NOSH_CSVWRITER_H
#define NOSH_CSVWRITER_H

#include <Teuchos_ParameterList.hpp>

namespace Nosh
{

class CsvWriter
{
public:

//! Default constructor.
  CsvWriter(const std::string &fileName,
            const std::string &delimeter = ","
          );

//! Destructor.
  virtual
  ~CsvWriter();

//! Const getter.

  void
  writeHeader(const Teuchos::ParameterList & pList) const;

  void
  writeRow(const Teuchos::ParameterList & pList) const;

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
} // namespace Nosh

#endif // NOSH_CSVWRITER_H
