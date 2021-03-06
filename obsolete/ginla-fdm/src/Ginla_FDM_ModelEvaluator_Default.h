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

#ifndef GINLA_FDM_MODELEVALUATOR_DEFAULT_H
#define GINLA_FDM_MODELEVALUATOR_DEFAULT_H

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <EpetraExt_ModelEvaluator.h>

#include "Ginla_FDM_Operator_Virtual.h"
#include "Komplex2_LinearProblem.h"
#include "Ginla_StateTranslator_Virtual.h"
#include "Ginla_CreateSavable_Virtual.h"
#include "Recti_Grid_General.h"

// forward declarations
namespace Ginla {
  namespace Operator {
    class Virtual;
  }
}

namespace Ginla {
namespace FDM {
namespace ModelEvaluator {

class Default:
    public EpetraExt::ModelEvaluator,
    public Ginla::StateTranslator::Virtual,
    public Ginla::CreateSavable::Virtual
{
public:

  //! Constructor without initial guess.
  Default ( const Teuchos::RCP<Ginla::FDM::Operator::Virtual> & glOperator,
            const Teuchos::RCP<Komplex2::LinearProblem>       & komplex,
            const Teuchos::RCP<Recti::Grid::General>          & grid,
            const Teuchos::ParameterList                      & params
          );

  //! Constructor with initial guess.
  Default ( const Teuchos::RCP<Ginla::FDM::Operator::Virtual> & glOperator,
            const Teuchos::RCP<Komplex2::LinearProblem>       & komplex,
            const Ginla::State::Updatable                     & state,
            const Teuchos::ParameterList                      & params
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

  virtual
  Teuchos::RCP<Ginla::State::Updatable>
  createState( const Epetra_Vector & x ) const;

  virtual
  Teuchos::RCP<Ginla::State::Virtual>
  createSavable( const Epetra_Vector & x ) const;

  virtual
  Teuchos::RCP<Epetra_Vector>
  createSystemVector( const Ginla::State::Updatable & state ) const;

  virtual
  void
  createSystemVector( const Ginla::State::Updatable & state,
                            Epetra_Vector           & x
                    ) const;

  Teuchos::RCP<const Epetra_Map>
  getMap() const;

public:
  // TODO move the following functions to private as soon as Ginla::ModelEvaluator::Bordered doesn't
  //      use them anymore
  void
  computeF ( const Epetra_Vector & x,
                   Epetra_Vector & FVec ) const;

  void
  computeDFDh0 ( const Epetra_Vector & x,
                       Epetra_Vector & FVec ) const;

  void
  computeJacobian ( const Epetra_Vector & x,
                    Epetra_Operator     & Jac ) const;

  Teuchos::RCP<Epetra_CrsMatrix>
  getJacobianNonConst();

public:
  // TODO remove the following functions as soon as Ginla::ModelEvaluator::Bordered doesn't
  //      use them anymore
  Teuchos::RCP<Ginla::FDM::Operator::Virtual>
  getOperator();

protected:
private:

   const Teuchos::RCP<Ginla::FDM::Operator::Virtual> glOperator_;
   const Teuchos::RCP<Komplex2::LinearProblem>       komplex_;
   const Teuchos::RCP<Epetra_Vector> x_;
   mutable bool firstTime_;

   int numParams_;

   Teuchos::RCP<Epetra_Map> p_map_;
   Teuchos::RCP<Epetra_Vector> p_init_;
   Teuchos::RCP<Teuchos::Array<std::string> > p_names_;

private:
    void
    setupParameters_( const Teuchos::ParameterList & params );

    Teuchos::RCP<Ginla::FDM::State>
    createFdmState_( const Epetra_Vector & x
                  ) const;

};
} // namespace ModelEvaluator
}
} // namespace Ginla

#endif // GINLA_FDM_MODELEVALUATOR_DEFAULT_H
