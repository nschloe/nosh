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

#ifndef VIO_IMAGE_WRITER_VTI_H
#define VIO_IMAGE_WRITER_VTI_H

#include "VIO_Image_Reader_Abstract.h"

namespace VIO
{
  namespace Image {
  namespace Reader
  {

class VTI:
    public VIO::Image::Reader::Abstract
{
public:
    VTI ( const std::string & filePath );

    virtual
    ~VTI ();

    virtual void
    read ( Teuchos::RCP<ComplexMultiVector>              & z,
           Teuchos::Array<int>                           & p,
           UIntTuple                                     & dims,
           Point                                         & origin,
           Point                                         & spacing,
           Teuchos::ParameterList                        & fieldData,
           const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
         ) const;

protected:
private:
};

  } // namespace Reader
 }
} // namespace VIO

#endif // VIO_IMAGE_WRITER_VTI_H
