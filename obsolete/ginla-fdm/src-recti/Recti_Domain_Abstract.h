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

#ifndef RECTI_DOMAIN_ABSTRACT_H
#define RECTI_DOMAIN_ABSTRACT_H

#include<Teuchos_Tuple.hpp>

typedef Teuchos::Tuple<double,3> Point;

namespace Recti
{
  namespace Domain
  {

class Abstract
{
public:
    //! Constructor.
    Abstract ( double tolerance = 1.0e-10 );

    //! Destructor.
    virtual
    ~Abstract();

    //! Returns a bounding box around the domain.
    virtual Teuchos::Tuple<double,4>
    getBoundingBox () const = 0;

    //! @param  x Point that is checked
    //! @return   Whether or not \c x sits in the domain.
    virtual bool
    isInDomain ( const Point & x ) const = 0;

protected:
    double tolerance_;
   
private:
};

  } // namespace Domain
} // namespace Recti

#endif // RECTI_DOMAIN_ABSTRACT_H
