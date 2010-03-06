/*
 * GridReader.h
 *
 *  Created on: Nov 27, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef GRIDREADER_H_
#define GRIDREADER_H_

#include "GridUniform.h"

#include <Tpetra_Vector.hpp>
#include <LOCA_Parameter_Vector.H>


namespace GridReader
{
void
read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
       const std::string                             & filePath,
       Teuchos::RCP<DoubleMultiVector>               & x,
       Teuchos::RCP<GridUniform>                     & grid,
       Teuchos::ParameterList                        & params
     );

void
read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
       const std::string                             & filePath,
       Teuchos::RCP<ComplexMultiVector>              & z,
       Teuchos::RCP<GridUniform>                     & grid,
       Teuchos::ParameterList                        & params
     );
};

#endif /* GRIDREADER_H_ */
