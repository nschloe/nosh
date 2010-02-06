/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

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

#include "DomainVirtual.h"
#include <Teuchos_Array.hpp>

//! Defer this class until the availabilty of GGL through Boost, providing polygon algorithms.
class DomainPolygon: public DomainVirtual
{
public:
    //! Constructor.
    DomainPolygon ( const Teuchos::Array<Teuchos::Tuple<double,2> > & polygonPoints,
                    const Teuchos::Tuple<double,2>                  & anchorPoint
                  );

    //! Destructor.
    virtual
    ~DomainPolygon();

    //! \brief  Classical point-in-polygon problem (http://en.wikipedia.org/wiki/Point_in_polygon).
    //! @param  x Point that is checked
    //! @return   Whether or not \c x sits in the domain.
    virtual bool
    isInDomain ( const Teuchos::Array<double> & x ) const;

protected:
private:
    Teuchos::Array<Teuchos::Tuple<double,2> > polygonPoints_;
};

#endif // DOMAINPOLYGON_H
