/*
 * GridReader.h
 *
 *  Created on: Nov 27, 2009
 *      Author: Nico Schlšmer
 */

#ifndef GRIDREADER_H_
#define GRIDREADER_H_

#include "GridVirtual.h"
#include "GridSquare.h"

class GridReader
{

public:
  GridReader();

  virtual
  ~GridReader();


  void
  read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
        const std::string                             & filePath,
        Teuchos::RCP<DoubleMultiVector>               & x,
        Teuchos::RCP<GridUniformVirtual>              & grid,
        Teuchos::ParameterList                        & params
      ) const;

protected:
private:

};

#endif /* GRIDREADER_H_ */
