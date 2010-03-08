/*
 * GridReader.h
 *
 *  Created on: Nov 27, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef GRIDREADER_H_
#define GRIDREADER_H_

#include "Recti_Grid_Uniform.h"

#include <Tpetra_Vector.hpp>
#include <LOCA_Parameter_Vector.H>

namespace Recti
{
  namespace Grid
  {
namespace Reader
{
void
read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
       const std::string                             & filePath,
       Teuchos::RCP<DoubleMultiVector>               & x,
       Teuchos::RCP<Uniform>                         & grid,
       Teuchos::ParameterList                        & params
     );

void
read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
       const std::string                             & filePath,
       Teuchos::RCP<ComplexMultiVector>              & z,
       Teuchos::RCP<Uniform>                         & grid,
       Teuchos::ParameterList                        & params
     );
};

  } // namespace Grid
} // namespace Reader

#endif /* GRIDREADER_H_ */
