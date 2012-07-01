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

#ifndef GINLA_CREATESAVABLE_VIRTUAL_H
#define GINLA_CREATESAVABLE_VIRTUAL_H

#include "Ginla_State_Virtual.h"
#include <Epetra_Vector.h>

namespace Ginla {
namespace CreateSavable {

class Virtual
{
public:
    Virtual();

    virtual
    ~Virtual();

    //! Translates a system vector into a state.
    virtual
    Teuchos::RCP<Ginla::State::Virtual>
    createSavable( const Epetra_Vector & x ) const = 0;

protected:
private:
};

} // namespace CreateSavable
} // namespace Ginla

#endif // GINLA_CREATESAVABLE_VIRTUAL_H
