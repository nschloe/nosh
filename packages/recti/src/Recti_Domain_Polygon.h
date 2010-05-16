/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

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

#ifndef DOMAINPOLYGON_H
#define DOMAINPOLYGON_H

#include "Recti_Domain_Abstract.h"
#include <Teuchos_Array.hpp>

namespace Recti
{
  namespace Domain
  {

//! Defer this class until the availabilty of GGL through Boost, providing polygon algorithms.
class Polygon:
   public Abstract
{
public:
    //! Constructor.
    Polygon ( const Teuchos::Array<DoubleTuple> & polygonPoints );

    //! Destructor.
    virtual
    ~Polygon();

    Teuchos::Tuple<double,4>
    getBoundingBox () const;

    //! \brief  Classical point-in-polygon problem (http://en.wikipedia.org/wiki/Point_in_polygon).
    //! @param  x Point that is checked
    //! @return   Whether or not \c x sits in the domain.
    virtual bool
    isInDomain ( const DoubleTuple & x ) const;

protected:
private:
    Teuchos::Array<DoubleTuple> polygonPoints_;
};

  } // namespace Domain
} // namespace Recti

#endif // DOMAINPOLYGON_H
