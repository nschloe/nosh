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

#ifndef VIO_READER_ABSTRACT_H
#define VIO_READER_ABSTRACT_H

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

namespace VIO 
{
  namespace Reader 
  {

class Abstract
{
public:
    Abstract ( const std::string & filePath );

    virtual
    ~Abstract ();

protected:
    int
    extractIntValue_ ( const vtkSmartPointer<vtkDataArray> & array
                     ) const;

    Teuchos::Array<int>
    extractIntArray_ ( const vtkSmartPointer<vtkDataArray> & array
                     ) const;

    double
    extractDoubleValue_ ( const vtkSmartPointer<vtkDataArray> & array
                        ) const;

    Teuchos::Array<double>
    extractDoubleArray_ ( const vtkSmartPointer<vtkDataArray> & array
                        ) const;

protected:

    Teuchos::ParameterList
    readFieldData_ ( const vtkSmartPointer<vtkDataObject>  & dataObject
                   ) const;

    std::string filePath_;
                  
private:
};

  } // namespace Reader
} // namespace VIO

#endif // VIO_READER_ABSTRACT_H
