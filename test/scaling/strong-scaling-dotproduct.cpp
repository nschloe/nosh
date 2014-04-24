// @HEADER
//
//    Scalability tests.
//    Copyright (C) 2012--2014  Nico Schl\"omer
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
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_VerboseObject.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

// =============================================================================
int myPow(int n, int k)
{
  int out = n;
  for (int kk = 0; kk < k; k++)
    out *= n;
  return out;
}
// =============================================================================
int main (int argc, char *argv[])
{
  // Initialize MPI
  Teuchos::GlobalMPISession (&argc, &argv, NULL);

  // Create output stream. (Handy for multicore output.)
  const RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  RCP<Epetra_MpiComm> eComm =
    rcp<Epetra_MpiComm> (new Epetra_MpiComm (MPI_COMM_WORLD));
#else
  RCP<Epetra_SerialComm>  eComm =
    rcp<Epetra_SerialComm> (new Epetra_SerialComm());
#endif

  bool success = true;
  try {
    // Create map.
    // Do strong scaling tests, so keep numGlobalElements independent of
    // the number of processes.
    const int maxSize = 9;
    for (int k = 0; k != maxSize+1; k++) {
      int numGlobalElements = myPow(10, k);
      // create map
      int indexBase = 0;
      RCP<Epetra_Map> map =
        rcp(new Epetra_Map (numGlobalElements, indexBase, *eComm));
      // create vectors
      RCP<Epetra_Vector> u = rcp(new Epetra_Vector(*map));
      u->Random();
      RCP<Epetra_Vector> v = rcp(new Epetra_Vector(*map));
      v->Random();

      // create timer label
      std::stringstream label;
      label << "Vector::Dot 10^" << k;

      RCP<Teuchos::Time> dotTime =
        Teuchos::TimeMonitor::getNewTimer(label.str());
      {
        Teuchos::TimeMonitor tm(*dotTime);
        double dot;
        TEUCHOS_ASSERT_EQUALITY(0, u->Dot(*v, &dot));
      }
    }
    // print timing data
    Teuchos::TimeMonitor::summarize();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
