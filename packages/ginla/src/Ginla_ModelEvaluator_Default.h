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

#ifndef GINLA_MODELEVALUATOR_GINLA_MODELEVALUATOR_DEFAULT_H
#define GINLA_MODELEVALUATOR_GINLA_MODELEVALUATOR_DEFAULT_H

#include <EpetraExt_ModelEvaluator.h>

#include "Ginla_Operator_Virtual.h"
#include "Ginla_Komplex.h"

#include <Thyra_ModelEvaluator.hpp>

namespace Ginla {

namespace ModelEvaluator {

class Default:
    public EpetraExt::ModelEvaluator
{
public:

  //! Constructor without initial guess.
  Default ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
            const Teuchos::RCP<Ginla::Komplex>           & komplex
          );

  //! Constructor with initial guess.
  Default ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
            const Teuchos::RCP<Ginla::Komplex>           & komplex,
            const Teuchos::RCP<const ComplexVector>      & psi
          );
  
  // Destructor
  virtual
  ~Default();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  virtual
  Teuchos::RCP<const Epetra_Map>
  get_x_map() const;

  virtual
  Teuchos::RCP<const Epetra_Map>
  get_f_map() const;

  virtual
  Teuchos::RCP<const Epetra_Vector>
  get_x_init () const;
  
  virtual
  Teuchos::RCP<Epetra_Operator> 
  create_W() const;

  virtual
  InArgs
  createInArgs() const;

  virtual
  OutArgs
  createOutArgs() const;

  virtual
  void
  evalModel( const InArgs  & inArgs,
             const OutArgs & outArgs ) const;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -              
  Teuchos::RCP<const Ginla::Komplex>
  getKomplex() const;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

protected:
private:
  
   const Teuchos::RCP<Ginla::Operator::Virtual> glOperator_;
   const Teuchos::RCP<Ginla::Komplex> komplex_;
//    const Teuchos::RCP<const Epetra_Vector> x_;
   const Teuchos::RCP<Epetra_Vector> x_;
   mutable bool firstTime_;
   
private:
    void
    computeF_ ( const Epetra_Vector & x,
                Epetra_Vector       & FVec ) const;
                
    void
    computeJacobian_ ( const Epetra_Vector & x,
                       Epetra_Operator     & Jac ) const;

};
} // namespace ModelEvaluator
} // namespace Ginla

#endif // GINLA_MODELEVALUATOR_GINLA_MODELEVALUATOR_DEFAULT_H
