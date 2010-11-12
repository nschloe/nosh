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

#ifndef GINLA_FDM_MODELEVALUATOR_BORDERED_H
#define GINLA_FDM_MODELEVALUATOR_BORDERED_H

#include <EpetraExt_ModelEvaluator.h>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Import.h>

#include "Ginla_StateTranslator_Virtual.h"
#include "Ginla_FDM_ModelEvaluator_Default.h"

// forward declarations
namespace Ginla {
  namespace Komplex {
    class LinearProblem;
  }
  namespace Operator {
    class Virtual;
  }
}

namespace Ginla {
namespace FDM {
namespace ModelEvaluator {

class Bordered:
    public EpetraExt::ModelEvaluator,
    public Ginla::StateTranslator::Virtual
{
public:

  //! Constructor without initial guess.
  Bordered ( const Teuchos::RCP<Ginla::FDM::Operator::Virtual> & glOperator,
             const Teuchos::RCP<Komplex2::LinearProblem>       & komplex,
             const Teuchos::RCP<Recti::Grid::General>          & grid,
             const Teuchos::ParameterList                      & params
           );

  //! Constructor with initial guess.
  Bordered ( const Teuchos::RCP<Ginla::FDM::Operator::Virtual> & glOperator,
             const Teuchos::RCP<Komplex2::LinearProblem>       & komplex,
             const Ginla::State::Updatable                     & state,
             const Teuchos::ParameterList                      & params
           );

  // Destructor
  virtual
  ~Bordered();

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
  Teuchos::RCP<Epetra_Vector>
  createSystemVector( const Ginla::State::Updatable & state ) const;

  virtual
  void
  createSystemVector( const Ginla::State::Updatable & state,
                            Epetra_Vector           & x
                    ) const;

protected:
private:

   const Teuchos::RCP<Ginla::FDM::ModelEvaluator::Default> modelEvaluatorDefault_;

   mutable bool firstTime_;

    Teuchos::RCP<const Epetra_BlockMap> regularMap_;
    Teuchos::RCP<Epetra_Map> extendedMap_;

    const Epetra_Import importFromExtendedMap_;
    const Epetra_Import importFromRegularMap_;

    const Teuchos::RCP<Epetra_CrsMatrix> jacobian_;
    const Teuchos::RCP<Epetra_Vector> x_;

private:

    void
    computeF_ ( const Epetra_Vector & x,
                      Epetra_Vector & FVec ) const;

    void
    computeDFDh0_ ( const Epetra_Vector & x,
                          Epetra_Vector & FVec ) const;

    void
    computeJacobian_ ( const Epetra_Vector & x,
                       Epetra_Operator     & Jac ) const;

    Teuchos::RCP<Epetra_Map>
    createExtendedMap_ ( const Epetra_BlockMap & realMap ) const;


    void
    fillBorderedMatrix_ ( const Teuchos::RCP<      Epetra_CrsMatrix> & extendedMatrix,
                          const Teuchos::RCP<const Epetra_CrsMatrix> & regularMatrix,
                          const Epetra_Vector                        & rightBorder,
                          // TODO Declare the following const as soon as Trilinos allows (ReplaceGlobalValues)
                          Epetra_Vector                              & lowerBorder,
                          double                                       d,
                          bool                                         firstTime
                        ) const;

    int
    PutRow_ ( const Teuchos::RCP<Epetra_CrsMatrix> A,
              const int                            globalRow,
              const int                            numIndices,
              double                             * values,
              int                                * indices,
              const bool                           firstTime
            ) const;
};

}
}
}

#endif // GINLA_FDM_MODELEVALUATOR_BORDERED_H
