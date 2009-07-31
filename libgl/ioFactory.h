#include "ioVirtual.h"

class IoFactory
{

  public:
     static IoVirtual* createFileIo( std::string fileName );

};