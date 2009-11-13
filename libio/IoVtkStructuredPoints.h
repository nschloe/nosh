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
		 ) const;

	virtual void
	write( const Tpetra::MultiVector<double,int> & x,
		   const int                               Nx,
		   const double                            h,
           const Teuchos::ParameterList          & problemParams
		 ) const;

protected:
private:

	void
	createOneProcMap( const Tpetra::Map<int> & sourceMap,
			                Teuchos::RCP<Tpetra::Map<int> > & oneProcMap );

	int exportProc_;
	Teuchos::RCP<Tpetra::Map<int> > sourceMap_;
	Teuchos::RCP<Tpetra::Export<int> >    oneProcExporter_;

};

#endif /* IOVTKSTRUCTUREDPOINTS_H_ */
