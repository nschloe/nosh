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

#ifndef GINLA_MODELEVALUATOR_DEFAULT_H
#define GINLA_MODELEVALUATOR_DEFAULT_H

#include <EpetraExt_ModelEvaluator.h>

#include "Ginla_Komplex.h"
#include "Ginla_State.h"

// forward declarations
namespace Ginla {
  namespace Operator {
    class Virtual;
  }
}

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
            const Teuchos::RCP<const Ginla::State>       & state
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
  Teuchos::RCP<const Epetra_Vector>
  get_p_init ( int l ) const;
  
  virtual
  Teuchos::RCP<const Epetra_Map>
  get_p_map(int l) const;
  
  virtual
  Teuchos::RCP<const Teuchos::Array<std::string> >
  get_p_names (int l) const;
  
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

  Teuchos::RCP<Ginla::State>
  createState( const Epetra_Vector & x ) const;

protected:
private:
  
   const Teuchos::RCP<Ginla::Operator::Virtual> glOperator_;
   const Teuchos::RCP<Ginla::Komplex> komplex_;
   const Teuchos::RCP<Epetra_Vector> x_;
   mutable bool firstTime_;
   
   const int numParameters_;
   
   Teuchos::RCP<Epetra_Map> p_map_;
   Teuchos::RCP<Epetra_Vector> p_init_;
   Teuchos::RCP<const Teuchos::Array<std::string> > p_names_;
   
private:
    void
    computeF_ ( const Epetra_Vector & x,
                Epetra_Vector       & FVec ) const;
                
    void
    computeJacobian_ ( const Epetra_Vector & x,
                       Epetra_Operator     & Jac ) const;
    
    Teuchos::RCP<Epetra_Vector>
    createX_(  const Ginla::State & state ) const;

};
} // namespace ModelEvaluator
} // namespace Ginla

#endif // GINLA_MODELEVALUATOR_DEFAULT_H
