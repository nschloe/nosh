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

#ifndef VTKREADER_H
#define VTKREADER_H

#include "AbstractImageReader.h"

#include <vtkStructuredPoints.h>

class VtkReader:
            public AbstractImageReader
{
public:
    VtkReader ( const std::string & filePath );

    virtual
    ~VtkReader ();

    virtual void
    read ( Teuchos::RCP<ComplexMultiVector>              & z,
           Teuchos::Array<int>                           & p,
           UIntTuple                                     & dims,
           DoubleTuple                                   & origin,
           DoubleTuple                                   & spacing,
           Teuchos::ParameterList                        & fieldData,
           const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
         ) const;

protected:
private:
};

#endif // VTKREADER_H
