/*
 * IoVtkStructuredPoints.cpp
 *
 *  Created on: Nov 10, 2009
 *      Author: Nico Schlšmer
 */

#include "IoVtkStructuredPoints.h"

//#include "vtkCell.h"
//#include "vtkCellData.h"
//#include "vtkStructuredGrid.h"
//#include "vtkStructuredGridWriter.h"
//#include "vtkImageData.h"
#include "vtkStructuredPoints.h"
//#include "vtkStructuredPointsReader.h"
#include "vtkStructuredPointsWriter.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"

// ============================================================================
IoVtkStructuredPoints::IoVtkStructuredPoints(std::string fname) :
  IoVirtual(fname), exportProc_(0), sourceMap_(Teuchos::null),
      oneProcExporter_(Teuchos::null)
{
}
// ============================================================================
IoVtkStructuredPoints::~IoVtkStructuredPoints()
{
}
// ============================================================================
void
IoVtkStructuredPoints::createOneProcMap(const Tpetra::Map<int> & sourceMap,
    Teuchos::RCP<Tpetra::Map<int> > & oneProcMap)
{
  int numGlobalElements = sourceMap.getGlobalNumElements();
  int numLocalElements;
  if (sourceMap.getComm()->getRank() == exportProc_)
    numLocalElements = numGlobalElements;
  else
    numLocalElements = 0;

  oneProcMap = Teuchos::rcp(new Tpetra::Map<int>(numGlobalElements,
      numLocalElements, sourceMap.getIndexBase(), sourceMap.getComm()));
}
// ============================================================================
void
IoVtkStructuredPoints::read(
    const Teuchos::RCP<const Teuchos::Comm<int> > &tComm, Teuchos::RCP<
        Tpetra::MultiVector<double, int> > &x,
    Teuchos::ParameterList &problemParams) const
{
  TEST_FOR_EXCEPTION( true,
      std::logic_error,
      "Not yet implemented." );
}
// ============================================================================
void
IoVtkStructuredPoints::write(const Tpetra::MultiVector<double, int> & x,
    const int Nx, const double h) const
{
  Teuchos::RCP<const Tpetra::Map<int> > sourceMapPtr = x.getMap();
  if (!x.getMap()->isSameAs(*sourceMap_))
    {
      // do a deep copy of the map
      //		Teuchos::rcp( x.getMap() ); // copy constructor is private!
      //		const Teuchos::RCP<const Tpetra::Map<int> > testMap( x.getMap() );
      sourceMap_ = x.getMap();
      // recreate the the exporter
      Teuchos::RCP<Tpetra::Map<int> > oneProcMap;
      createOneProcMap(*sourceMap_, oneProcMap);
      oneProcExporter_ = Teuchos::rcp(new Tpetra::Export<int>(*sourceMap_,
          *oneProcMap));
    }

  // export all the data to one processor
  const Tpetra::MultiVector<double, int> xOneProc;
  x.doExport(xOneProc, oneProcExporter_);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // All the date now sits on exportProc_; write the file from here.
  if (sourceMap_->getComm()->getRank() == exportProc_){
      vtkStructuredPoints * spData = vtkStructuredPoints::New();

      spData->SetDimensions(Nx + 1, Nx + 1, 1);
      spData->SetOrigin(0, 0, 0);
      spData->SetSpacing(h, h, 0);

      Teuchos::ArrayRCP<const double> xOneProcView =
          xOneProc.getVector(0)->get1dView();

      vtkDoubleArray * scalars = vtkDoubleArray::New();
      int k = 0;
      for (int z = 0; z < 1; z++)
        for (int y = 0; y < Nx + 1; y++)
          for (int x = 0; x < Nx + 1; x++)
            scalars->InsertNextValue(xOneProcView[k++]);

      spData->GetPointData()->SetScalars(scalars);

      // write the stuff
      vtkStructuredPointsWriter* writer = vtkStructuredPointsWriter::New();
      writer->SetInput(spData);
      writer->SetHeader("myheader");
      writer->SetScalarsName("myscalarsname");
      writer->SetFileName(fileName_.c_str());
      writer->SetFileTypeToASCII();
      //    writer->SetFileTypeToBinary();
      writer->Write();

      // clean up
      spData->Delete();
      scalars->Delete();
      writer->Delete();
    }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
// ============================================================================
void
IoVtkStructuredPoints::write(const Tpetra::MultiVector<double, int> & x,
    const int Nx, const double h, const Teuchos::ParameterList & problemParams) const
{
  write(x, Nx, h);
}
// =============================================================================
