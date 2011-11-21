// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2010, 2011  Nico Schl\"omer
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
#ifndef GINLA_MODELEVALUATOR_H
#define GINLA_MODELEVALUATOR_H
// -----------------------------------------------------------------------------
// includes
// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include "Ginla_config.h"

#include <EpetraExt_ModelEvaluator.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif

#include "Ginla_JacobianOperator.hpp"
// -----------------------------------------------------------------------------
// forward declarations
namespace Ginla
{
  class MagneticVectorPotential;
  class State;
  class StkMesh;
}

class Epetra_CrsGraph;
// -----------------------------------------------------------------------------
namespace Ginla {

class ModelEvaluator: public EpetraExt::ModelEvaluator
{

public:
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //! Constructor without initial guess.
  ModelEvaluator ( const Teuchos::RCP<Ginla::StkMesh>      & mesh,
                   const Teuchos::ParameterList                       & params,
                   const Teuchos::RCP<const Epetra_Vector>            & thickness,
                   const Teuchos::RCP<Ginla::MagneticVectorPotential> & mvp,
                   const Teuchos::RCP<Ginla::State>        & initialState
                 );

  // Destructor
  virtual
  ~ModelEvaluator();

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
  evalModel( const InArgs  & inArgs,
             const OutArgs & outArgs ) const;

  virtual
  Teuchos::RCP<Ginla::State>
  createSavable( const Epetra_Vector & x ) const ;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

public:

    virtual
    Teuchos::RCP<LOCA::ParameterVector>
    getParameters() const;

private:
  void
  computeF_ ( const Epetra_Vector                             & x,
              const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams,
              const double                                      temperature,
                    Epetra_Vector                             & FVec
            ) const;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
protected:
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
private:
//    const Teuchos::RCP<const Epetra_Comm> & eComm_;
//    Teuchos::RCP<const Epetra_Map>          eMap_;
//
//    const Teuchos::RCP<const TComm> & tComm_;
//    Teuchos::RCP<const TMap>          tMap_;

   const Teuchos::RCP<Ginla::StkMesh> mesh_;

   const Teuchos::RCP<const Epetra_Vector> thickness_;

   const Teuchos::RCP<Epetra_Vector> x_;

   int numParams_;

   Teuchos::RCP<Epetra_Map> p_map_;
   Teuchos::RCP<Epetra_Vector> p_init_;
   Teuchos::RCP<Teuchos::Array<std::string> > p_names_;

   // for get_parameters()
   Teuchos::RCP<Epetra_Vector> p_current_;

   const Teuchos::RCP<Ginla::MagneticVectorPotential> mvp_;

   const Teuchos::RCP<Ginla::KeoFactory> keoFactory_;

#ifdef GINLA_TEUCHOS_TIME_MONITOR
   const Teuchos::RCP<Teuchos::Time> evalModelTime_;
   const Teuchos::RCP<Teuchos::Time> computeFTime_;
   const Teuchos::RCP<Teuchos::Time> fillJacobianTime_;
   const Teuchos::RCP<Teuchos::Time> fillPreconditionerTime_;
#endif

  Teuchos::RCP<Teuchos::FancyOStream> out_;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
private:
    void
    setupParameters_( const Teuchos::ParameterList & params );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
};

}

#endif // GINLA_MODELEVALUATOR_H
