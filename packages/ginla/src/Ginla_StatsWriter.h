/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef GL_STATSWRITER_H
#define GL_STATSWRITER_H

#include <Teuchos_ParameterList.hpp>

namespace Ginla {

class StatsWriter
{
public:

  //! Default constructor.
  StatsWriter( std::string & fileName );
  
  //! Destructor.
  virtual
  ~StatsWriter();
  
  //! Non-const getter.
  Teuchos::RCP<Teuchos::ParameterList>
  getList();
  
  void
  setList( const Teuchos::RCP<Teuchos::ParameterList> statisticsList );
  
  void
  print();
  
protected:
private:
  
  //! File stream for the statistics.
  std::ofstream fileStream_;
  
  //! Stores the scalar statisticsList values
  Teuchos::RCP<Teuchos::ParameterList> statisticsList_;
  
  //! Whether or not to print the header in the stats file.
  bool printHeader_;
  
};

}

#endif // GL_STATSWRITER_H
