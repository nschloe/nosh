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

#ifndef DOMAINRECTANGLE_H
#define DOMAINRECTANGLE_H

#include "DomainVirtual.h"
#include <boost-1_35/boost/concept_check.hpp>

class DomainRectangle:
            public DomainVirtual
{
public:
    //! Constructor.
    DomainRectangle ( const double a,
                      const double b );

    //! Destructor.
    virtual
    ~DomainRectangle();

    //! Returns a bounding box around the domain.
    virtual Teuchos::Tuple<double,4>
    getBoundingBox () const;

    //! @param  x Point that is checked
    //! @return   Whether or not \c x sits in the domain.
    virtual bool
    isInDomain ( const Teuchos::Tuple<double,2> & x ) const;

protected:
private:
    double a_;
    double b_;
};

#endif // DOMAINRECTANGLE_H
