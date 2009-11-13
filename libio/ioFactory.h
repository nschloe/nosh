#ifndef IOFACTORY_H
#define IOFACTORY_H

#include "ioVirtual.h"

class IoFactory
{

public:
  //! Returns an pointer to an object of type IoVirtual for a given file
  //! fileName.
  //! Which of the implementations of IoVirtual is chosen is determined
  //! according to the suffix if fileName.
  static IoVirtual*
  createFileIo(std::string fileName);

};
#endif // IOFACTORY_H
