/*
 * GridReader.h
 *
 *  Created on: Nov 27, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef RECTI_GRID_READER_H_
#define RECTI_GRID_READER_H_

#include "Recti_Grid_Uniform.h"
namespace Ginla {
  class State;
};

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
       Teuchos::RCP<Ginla::State>                    & state,
       Teuchos::RCP<Recti::Grid::Uniform>            & grid,
       Teuchos::ParameterList                        & params
     );
};

  } // namespace Grid
} // namespace Reader

#endif /* RECTI_GRID_READER_H_ */
