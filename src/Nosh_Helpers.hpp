// @HEADER
//
//    Nosh helper functions.
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

#ifndef GL_HELPERS_H
#define GL_HELPERS_H

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <LOCA_Parameter_Vector.H>
#include <Teuchos_ParameterList.hpp>

namespace Nosh {

namespace Helpers {
//! Method that reads a given Teuchos::ParameterList and puts all \c double
//! entries into a LOCA::ParameterVector.
Teuchos::RCP<LOCA::ParameterVector>
teuchosParameterList2locaParameterVector( const Teuchos::ParameterList &p
                                          );

Teuchos::RCP<Teuchos::ParameterList>
locaParameterVector2teuchosParameterList( const LOCA::ParameterVector &pL );

//! Merges two \c LOCA::ParameterLists into one, checking for discrepancies
//! in the entries.
Teuchos::RCP<LOCA::ParameterVector>
mergeLocaParameterVectors( const LOCA::ParameterVector &p0,
                           const LOCA::ParameterVector &p1
                           );

void
appendToTeuchosParameterList( Teuchos::ParameterList &p,
                              const LOCA::ParameterVector &pL,
                              const std::string &labelPrepend = ""
                              );

bool
locaParameterVectorsEqual( const Teuchos::RCP<const LOCA::ParameterVector> &a,
                           const Teuchos::RCP<const LOCA::ParameterVector> &b
                           );

unsigned int
numDigits( const int i );

}

}

#endif // GL_HELPERS_H
