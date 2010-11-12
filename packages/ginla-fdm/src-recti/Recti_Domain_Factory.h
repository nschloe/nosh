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

#ifndef RECTI_DOMAIN_FACTORY_H
#define RECTI_DOMAIN_FACTORY_H

#include "Recti_Domain_Abstract.h"

#include <Teuchos_ParameterList.hpp>

namespace Recti {

namespace Domain {

class Factory
{
  public:
    //! constructor with paramlist
    Factory( const Teuchos::ParameterList & pList );
    
    //! Destructor.
    ~Factory();
    
    //! Construct the domain.
    Teuchos::RCP<Recti::Domain::Abstract>
    build() const;
    
  protected:
  private:
    Teuchos::RCP<Recti::Domain::Abstract>
    buildCircle() const;
    
    Teuchos::RCP<Recti::Domain::Abstract>
    buildEllipse() const;
    
    Teuchos::RCP<Recti::Domain::Abstract>
    buildPolygon() const;
    
    Teuchos::RCP<Recti::Domain::Abstract>
    buildRectangle() const;
    
    Teuchos::RCP<Recti::Domain::Abstract>
    buildSquare() const;
    
  private:
    Teuchos::ParameterList pList_;
};

// nonmember functions
Teuchos::RCP<Recti::Domain::Abstract>
buildDomain(  const Teuchos::ParameterList & pList );

}

}

#endif // RECTI_DOMAIN_FACTORY_H
