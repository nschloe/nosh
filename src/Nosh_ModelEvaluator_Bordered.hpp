// @HEADER
//
//    Nosh bordered model evaluator.
//    Copyright (C) 2012  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
#ifndef NOSH_MODELEVALUATOR_BORDERED_H
#define NOSH_MODELEVALUATOR_BORDERED_H
// -----------------------------------------------------------------------------
// includes
#include <Teuchos_RCP.hpp>

#include "Nosh_ModelEvaluator_Virtual.hpp"
// -----------------------------------------------------------------------------
namespace Nosh {
namespace ModelEvaluator {
class Bordered : public Nosh::ModelEvaluator::Virtual
{

public:
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//! Constructor without initial guess.
Bordered (
  const Teuchos::RCP<const Nosh::ModelEvaluator::Virtual> & modelEval,
  const Teuchos::RCP<const Epetra_Vector> & initialBordering,
  const double lambdaInit
  );

// Destructor
virtual
~Bordered();

virtual
Teuchos::RCP<const Epetra_Map>
get_x_map() const;

virtual
Teuchos::RCP<const Epetra_Map>
get_f_map() const;

virtual
Teuchos::RCP<const Epetra_Vector>
get_x_init() const;

virtual
Teuchos::RCP<const Epetra_Vector>
get_p_init( int l ) const;

virtual
Teuchos::RCP<const Epetra_Map>
get_p_map( int l ) const;

virtual
Teuchos::RCP<const Teuchos::Array<std::string> >
get_p_names( int l ) const;

virtual
Teuchos::RCP<Epetra_Operator>
create_W() const;

virtual
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
create_WPrec() const;

virtual
InArgs
createInArgs() const;

virtual
OutArgs
createOutArgs() const;

virtual
void
evalModel( const InArgs &inArgs,
           const OutArgs &outArgs ) const;

public:
virtual
Teuchos::RCP<const Epetra_Vector>
get_p_latest() const;

virtual
double
innerProduct(const Epetra_Vector &phi,
             const Epetra_Vector &psi
             ) const;

virtual
double
normalizedScaledL2Norm(const Epetra_Vector &psi) const;

virtual
double
gibbsEnergy(const Epetra_Vector &psi) const;

virtual
const Teuchos::RCP<const Nosh::StkMesh>
getMesh() const;

protected:
private:
const Teuchos::RCP<const Nosh::ModelEvaluator::Virtual> innerModelEval_;
const Teuchos::RCP<const Epetra_Vector> initialBordering_;
const double lambdaInit_;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
};
} // namespace Modelevaluator
} // namespace Nosh

#endif // NOSH_MODELEVALUATOR_BORDERED_H
