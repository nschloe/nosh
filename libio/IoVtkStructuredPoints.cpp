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
#include "vtkDataSetWriter.h"
//#include "vtkStructuredPointsReader.h"
#include "vtkStructuredPointsWriter.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

#include "vtkImageStencil.h"
#include "vtkXMLImageDataWriter.h"
#include "vtkImageAppendComponents.h"

// ============================================================================
IoVtkStructuredPoints::IoVtkStructuredPoints(std::string fname) :
  IoVirtual(fname), ioProc_(0), sourceMap_(Teuchos::null),
      oneProcImporter_(Teuchos::null)
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
  if (sourceMap.getComm()->getRank() == ioProc_)
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
    const int Nx, const double h)
{
  if ( !sourceMap_.is_valid_ptr() ||
        sourceMap_.is_null() ||
       !x.getMap()->isSameAs(*sourceMap_) )
    {
      sourceMap_ = x.getMap();
      // recreate the the exporter
      createOneProcMap(*sourceMap_, oneProcMap_);
      oneProcImporter_ = Teuchos::rcp(new Tpetra::Import<int>(sourceMap_,
          oneProcMap_));
    }

  // Export all the data to one processor.
  // Can't make it const as vtkDoubleArray::SetArray below needs non-const
  // access.
  // TODO See if this can be made const.
  unsigned int numVectors = x.getNumVectors();
  Teuchos::RCP<Tpetra::MultiVector<double, int> > xOneProc =
      Teuchos::rcp( new Tpetra::MultiVector<double,int>(oneProcMap_,numVectors) );
  xOneProc->doImport( x, *oneProcImporter_, Tpetra::INSERT );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // All the date now sits on ioProc_; write the file from here.
  if (sourceMap_->getComm()->getRank() == ioProc_){

    vtkImageAppendComponents * imageAppend = vtkImageAppendComponents::New();

    // TODO Replace by StructuredPoitnsWriter?
    // TODO Replace this function by a general data set preparator, and have one writer for them all?
    vtkSmartPointer<vtkDataSetWriter> writer = vtkDataSetWriter::New();
    writer->SetFileName(fileName_.c_str());
    writer->SetHeader("myheader");
    writer->SetFileTypeToASCII();
//    writer->SetFileTypeToBinary();
    writer->SetScalarsName("psi");

    for (unsigned int k=0; k<numVectors; k++ ) {
      vtkSmartPointer<vtkStructuredPoints> spData = vtkStructuredPoints::New();
      spData->SetDimensions(Nx + 1, Nx + 1, 1);
      spData->SetOrigin(0, 0, 0);
      spData->SetSpacing(h, h, 0);

      // NonConst view necessary here as vtkDoubleArray::SetArray needs
      // a nonconst array (for whatever reason).
      // TODO See if this can be made const.
      Teuchos::ArrayRCP<double> xOneProcView =
          xOneProc->getVectorNonConst(k)->get1dViewNonConst();

      // TODO Replace this by a length specification of xOneProcView
      int size = xOneProc->getGlobalLength();
      // TODO See if this can be made const double.
      double* array = xOneProcView.getRawPtr();
      int save = 1; // don't have VTK delete the array at cleanup
      vtkSmartPointer<vtkDoubleArray> scalars = vtkDoubleArray::New();
      scalars->SetArray( array, size, save );
      spData->GetPointData()->SetScalars(scalars);

      imageAppend->AddInput(spData);
    }

    writer->SetInputConnection(imageAppend->GetOutputPort());
    writer->Update(); // really necessary?
    writer->Write();
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
// ============================================================================
void
IoVtkStructuredPoints::write(const Tpetra::MultiVector<double, int> & x,
    const int Nx, const double h, const Teuchos::ParameterList & problemParams)
{
  write(x, Nx, h);
}
// =============================================================================
