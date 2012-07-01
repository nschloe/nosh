/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef VTKWRITER_H
#define VTKWRITER_H

#include "VIO_Image_Writer_Abstract.h"

namespace VIO
{
  namespace Image {
  namespace Writer
  {

class legacyVTK:
            public Abstract
{
public:

    //! Default constructor.
    legacyVTK ( const std::string & filePath );

    //! Destructor
    virtual ~legacyVTK();

    virtual void
    write () const;
    
    void
    setFormatBinary();
    
    void
    setFormatAscii();

protected:
private:
  //! whether to write binary of ASCII format VTK
  bool isFormatBinary_;
};

  } // namespace Writer
 }
} // namespace VIO

#endif // VTKWRITER_H
