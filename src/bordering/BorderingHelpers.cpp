// @HEADER
//
//    Nosh Helper functions.
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

#include "BorderingHelpers.hpp"

#include <vector>

#include <Teuchos::Comm<int>.h>
#include <Epetra_Import.h>

namespace nosh
{
// ============================================================================
std::shared_ptr<const Tpetra::Map<int,int>>
BorderingHelpers::
extendMapBy1(const Epetra_BlockMap & map)
{
  const Teuchos::Comm<int> & comm = map.Comm();
  // Create a new map that hosts one more entry.
  const int numGlobalElements = map.NumGlobalElements() + 1;
  const int numMyElements = map.NumMyElements();
  int * myGlobalElements = map.MyGlobalElements();
  // The following if-else construction just makes sure that
  // the Tpetra::Map<int,int> constructor is called with an extended
  // map on proc 0, and with the regular old stuff on all
  // other procs.
  std::shared_ptr<Tpetra::Map<int,int>> extendedMap;
  if (comm.MyPID() == 0) {
    // Copy over the global indices.
    std::vector<int> a(numMyElements+1);
    for (int k = 0; k < numMyElements; k++)
      a[k] = myGlobalElements[k];
    // Append one more.
    a[numMyElements] = map.NumGlobalElements();

    extendedMap = std::make_shared<Tpetra::Map<int,int>>(
        numGlobalElements,
        numMyElements+1,
        &a[0],
        map.IndexBase(),
        comm
        );
  } else {
    extendedMap = std::make_shared<Tpetra::Map<int,int>>(
        numGlobalElements,
        numMyElements,
        myGlobalElements,
        map.IndexBase(),
        comm
        );
  }

  return extendedMap;
}
// ============================================================================
void
BorderingHelpers::
merge(const Tpetra::MultiVector<double,int,int> & x,
      const double * lambda,
      Tpetra::MultiVector<double,int,int> & out
    )
{
#ifndef NDEBUG
  // Check if the maps are matching.
  std::shared_ptr<const Tpetra::Map<int,int>> extendedMap =
    nosh::BorderingHelpers::extendMapBy1(x.getMap());
  TEUCHOS_ASSERT(out.getMap().SameAs(*extendedMap));
#endif

  Epetra_Import importer(out.getMap(), x.getMap());

  TEUCHOS_ASSERT_EQUALITY(0, out.Import(x, importer, Insert));

  // Set last entry on proc 0.
  if (x.getMap().Comm().MyPID() == 0) {
    const int numMyElems = x.getMap().NumMyElements();
    for (int k = 0; k < x.NumVectors(); k++)
      (*out(k))[numMyElems] = lambda[k];
  }

  return;
}
// ============================================================================
void
BorderingHelpers::
dissect(const Tpetra::MultiVector<double,int,int> & x,
        Tpetra::MultiVector<double,int,int> & xSmall,
        double * lambda
       )
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(x.NumVectors(), xSmall.NumVectors());
  // Make sure the maps are matching.
  std::shared_ptr<const Tpetra::Map<int,int>> extendedMap =
    nosh::BorderingHelpers::extendMapBy1(xSmall.getMap());
  TEUCHOS_ASSERT(x.getMap().SameAs(*extendedMap));
#endif

  Epetra_Import importer(xSmall.getMap(), x.getMap());

  // Strip off the phase constraint variable.
  xSmall.Import(x, importer, Insert);

  // TODO Check if we need lambda on all procs.
  if (x.getMap().Comm().MyPID() == 0) {
    const int n = x.MyLength();
    for (int k = 0; k < x.NumVectors(); k++)
      lambda[k] = (*(x(k)))[n - 1];
  }

  return;
}
// ============================================================================
} // namespace nosh
