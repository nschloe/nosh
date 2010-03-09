/*
 * AbstractStateWriter.h
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef ABSTRACTSTATEWRITER_H_
#define ABSTRACTSTATEWRITER_H_

#include <Epetra_Vector.h>

//! Purely abstract base class that provides functionality for writing
//! states (e.g., solutions) of a system conneted with a geometry.
class AbstractStateWriter
{
public:
  //! Constructor
  AbstractStateWriter();

  // Destructor.
  virtual
  ~AbstractStateWriter();

  //! Method that writes the solution into the file \c filename.
  virtual void
  writeSolutionToFile( const Epetra_Vector & solution,
                       const std::string   & fileName ) const = 0;

  //! Method that writes an abstract state into the file \c filename.
  virtual void
  writeAbstractStateToFile( const Epetra_Vector & solution,
                            const std::string   & fileName ) const = 0;

protected:
private:
};

#endif /* ABSTRACTSTATEWRITER_H_ */
