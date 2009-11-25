/*
 * IoVtkStructuredPoints.h
 *
 *  Created on: Nov 10, 2009
 *      Author: Nico Schlšmer
 */

#ifndef IOVTKSTRUCTUREDPOINTS_H_
#define IOVTKSTRUCTUREDPOINTS_H_

#include "ioVirtual.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Export.hpp>

#include <vtkLookupTable.h>
#include <vtkSmartPointer.h>

class IoVtkStructuredPoints: public IoVirtual
{
public:

	IoVtkStructuredPoints( std::string fname );

	virtual ~IoVtkStructuredPoints();

	virtual void
	read( const Teuchos::RCP<const Teuchos::Comm<int> >        &tComm,
	            Teuchos::RCP<Tpetra::MultiVector<double,int> > &x,
	            Teuchos::ParameterList                         &problemParams
	    ) const;


	virtual void
	write( const Tpetra::MultiVector<double,int> & x,
		   const int                               Nx,
		   const double                            h
		 );

	virtual void
	write( const Tpetra::MultiVector<double,int> & x,
		   const int                               Nx,
		   const double                            h,
           const Teuchos::ParameterList          & problemParams
		 );

protected:
private:

	void
	createOneProcMap( const Tpetra::Map<int>          & sourceMap,
                          Teuchos::RCP<Tpetra::Map<int> > & oneProcMap );

	void
	constructLookupTable(vtkSmartPointer<vtkLookupTable> LookupTable) const;

	int ioProc_;
	Teuchos::RCP<const Tpetra::Map<int> > sourceMap_;
	Teuchos::RCP<Tpetra::Map<int> >       oneProcMap_;
	Teuchos::RCP<Tpetra::Import<int> > oneProcImporter_;

};

#endif /* IOVTKSTRUCTUREDPOINTS_H_ */
