// @HEADER
//
//    Nosh helper functions for bordered systems.
//    Copyright (C) 2010--2012  Nico Schlömer
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

#ifndef NOSH_BORDERINGHELPERS_H
#define NOSH_BORDERINGHELPERS_H

#include <Teuchos_RCP.hpp>
#include <Tpetra::Map<int,int>.h>
#include <Tpetra_Vector.hpp>

namespace Nosh
{
namespace BorderingHelpers
{
std::shared_ptr<const Tpetra::Map<int,int>>
extendMapBy1(const Epetra_BlockMap & map);

void
merge(const Tpetra::MultiVector<double,int,int> & x,
      const double * lambda,
      Tpetra::MultiVector<double,int,int> & out
    );

void
dissect(const Tpetra::MultiVector<double,int,int> & x,
        Tpetra::MultiVector<double,int,int> & xSmall,
        double * lambda
      );
} // namespace BorderingHelpers
} // namespace Nosh
#endif // NOCH_BORDERINGHELPERS_H
