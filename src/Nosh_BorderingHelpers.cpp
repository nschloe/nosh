// @HEADER
//
//    Nosh Helper functions.
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

#include "Nosh_BorderingHelpers.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Import.h>

namespace Nosh {
// ============================================================================
Teuchos::RCP<const Epetra_Map>
BorderingHelpers::
extendMapBy1(const Epetra_BlockMap & map)
{
  const Epetra_Comm & comm = map.Comm();
  // Create a new map that hosts one more entry.
  const int numGlobalElements = map.NumGlobalElements()
                              + 1;
  int numMyElements = map.NumMyElements();
  int * myGlobalElements = map.MyGlobalElements();
  // The following if-else construction just makes sure that
  // the Epetra_Map constructor is called with an extended
  // map on proc 0, and with the regular old stuff on all
  // other procs.
  Teuchos::RCP<Epetra_Map> extendedMap;
  if (comm.MyPID() == 0)
  {
    // Copy over the global indices.
    int a[numMyElements+1];
    for (int k=0; k<numMyElements; k++)
      a[k] = myGlobalElements[k];
    // Append one more.
    a[numMyElements] = map.NumGlobalElements();

    extendedMap =
      Teuchos::rcp(new Epetra_Map(numGlobalElements,
                                  numMyElements+1,
                                  a,
                                  map.IndexBase(),
                                  comm));
  }
  else
  {
    extendedMap =
      Teuchos::rcp(new Epetra_Map(numGlobalElements,
                                  numMyElements,
                                  myGlobalElements,
                                  map.IndexBase(),
                                  comm));
  }

  return extendedMap;
}
// ============================================================================
void
BorderingHelpers::
merge(const Epetra_MultiVector & x,
      const double * lambda,
      Epetra_MultiVector & out
      )
{
#if _DEBUG_
  TEUCHOS_ASSERT_EQUALITY(x.Map().NumGlobalEntries() + 1,
                          out.Map().NumGlobalEntries());
#endif
  Teuchos::RCP<const Epetra_Map> extendedMap =
    Nosh::BorderingHelpers::extendMapBy1(x.Map());

  Epetra_Import importer(*extendedMap, x.Map());

#if _DEBUG_
  out.Map.IsSameAs(*extendedMap);
#endif
  TEUCHOS_ASSERT_EQUALITY(0, out.Import(x, importer, Insert));

  // Set last entry on proc 0.
  if (x.Map().Comm().MyPID() == 0)
  {
    const int numMyElems = x.Map().NumMyElements();
    for (int k=0; k<x.NumVectors(); k++)
      (*out(k))[numMyElems] = lambda[k];
  }

  return;
}
// ============================================================================
void
BorderingHelpers::
dissect(const Epetra_MultiVector & x,
        Epetra_MultiVector & xSmall,
        double * lambda
        )
{
#if _DEBUG_
  TEUCHOS_ASSERT_EQUALITY(x.NumVectors(), xSmall.NumVectors());
#endif
  Teuchos::RCP<const Epetra_Map> extendedMap =
    Nosh::BorderingHelpers::extendMapBy1(xSmall.Map());
  Epetra_Import importer(xSmall.Map(),
                         *extendedMap);
#if _DEBUG_
  TEUCHOS_ASSERT(x.Map.IsSameAs(*extendedMap));
#endif

  // Strip off the phase constraint variable.
  xSmall.Import(x, importer, Insert);

  // TODO local access to global element ugh!
  for (int k=0; k<x.NumVectors(); k++)
    lambda[k] = (*(x(k)))[x.GlobalLength() - 1];

  return;
}
// ============================================================================
} // namespace Nosh
