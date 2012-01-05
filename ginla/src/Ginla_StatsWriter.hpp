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
#ifndef GINLA_STATSWRITER_H
#define GINLA_STATSWRITER_H

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

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
  
  //! Const getter.
  Teuchos::RCP<const Teuchos::ParameterList>
  getList();
  
  //! Non-const getter.
  Teuchos::RCP<Teuchos::ParameterList>
  getListNonConst();
  
  void
  setList( const Teuchos::ParameterList & statisticsList );
  
  void
  setList( const Teuchos::RCP<Teuchos::ParameterList> & statisticsList );
  
  void
  print();
  
protected:
private:
  
  //! File stream for the statistics.
  std::ofstream fileStream_;
  
  //! Stores the scalar statisticsList values
  const Teuchos::RCP<Teuchos::ParameterList> statisticsList_;
  
  //! Whether or not to print the header in the stats file.
  bool printHeader_;
  
};
} // namespace Ginla

#endif // GINLA_STATSWRITER_H
