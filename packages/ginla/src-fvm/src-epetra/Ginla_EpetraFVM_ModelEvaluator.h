/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schll\"omer

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

#ifndef GINLA_EPETRAFVM_EVALUATOR_H
#define GINLA_EPETRAFVM_EVALUATOR_H
// -----------------------------------------------------------------------------
// includes
#include <EpetraExt_ModelEvaluator.h>
#include <Epetra_Vector.h>


#include "VIO_EpetraMesh_Mesh.h"
#include "Ginla_ParameterHost_Virtual.h"
#include "Ginla_EpetraFVM_JacobianOperator.h"
#include "Ginla_CreateSavable_Virtual.h"
// -----------------------------------------------------------------------------
//typedef Tpetra::CrsGraph<ORD> TCrsGraph;
//typedef Tpetra::CrsMatrix<double_complex,ORD> ComplexMatrix;
// -----------------------------------------------------------------------------
// forward declarations
namespace Ginla
{
  namespace MagneticVectorPotential
  {
    class Virtual;
  }
  namespace EpetraFVM {
    class State;
  }
}

class Epetra_CrsGraph;
// -----------------------------------------------------------------------------
namespace Ginla {
namespace EpetraFVM {

class ModelEvaluator:
    public EpetraExt::ModelEvaluator,
    public Ginla::ParameterHost::Virtual,
    public Ginla::CreateSavable::Virtual
{

public:
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //! Constructor without initial guess.
  ModelEvaluator ( const Teuchos::RCP<VIO::EpetraMesh::Mesh>                   & mesh,
                   const Teuchos::ParameterList                                & params,
                   const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp,
                   const Teuchos::RCP<Ginla::EpetraFVM::State>                 & initialState
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
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

public:

    virtual
    Teuchos::RCP<LOCA::ParameterVector>
    getParameters() const;

private:
  void
  computeF_ ( const Epetra_Vector            & x,
              const double                     lambda,
              const Teuchos::Tuple<double,3> & scaling,
              const double                     temperature,
                    Epetra_Vector            & FVec
            ) const;

//  void
//  computeDFDp_ ( const Epetra_Vector            & x,
//                 const double                     mu,
//                 const Teuchos::Tuple<double,3> & scaling,
//                       Epetra_Vector            & FVec
//               ) const;

  void
  computeJacobian_ ( const Teuchos::RCP<const Epetra_Vector> & x,
                     const double                              lambda,
                     const Teuchos::Tuple<double,3>          & scaling,
                     const double                              temperature,
                           Epetra_Operator                   & Jac
                   ) const;

  void
  computePreconditioner_ ( const double                     mu,
                           const Teuchos::Tuple<double,3> & scaling,
                           Epetra_Operator                & Prec
                         ) const;

  virtual
  Teuchos::RCP<Ginla::State::Virtual>
  createSavable( const Epetra_Vector & x ) const ;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
protected:
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
private:
//    const Teuchos::RCP<const Epetra_Comm> & eComm_;
//    Teuchos::RCP<const Epetra_Map>          eMap_;
//
//    const Teuchos::RCP<const TComm> & tComm_;
//    Teuchos::RCP<const TMap>          tMap_;

   const Teuchos::RCP<VIO::EpetraMesh::Mesh> mesh_;

   Teuchos::RCP<Epetra_Vector> x_;
   mutable bool firstTime_;

   int numParams_;

   Teuchos::RCP<Epetra_Map> p_map_;
   Teuchos::RCP<Epetra_Vector> p_init_;
   Teuchos::RCP<Teuchos::Array<std::string> > p_names_;

   // for get_parameters()
   Teuchos::RCP<Epetra_Vector> p_current_;

   const Teuchos::RCP<Ginla::EpetraFVM::KineticEnergyOperator> keo_;
   const Teuchos::RCP<Ginla::EpetraFVM::JacobianOperator> jacobianOperator_;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
private:
    void
    setupParameters_( const Teuchos::ParameterList & params );

    // this value is really just for getParameters()
    mutable double mu_;
    mutable double scaling_;
    mutable Teuchos::Tuple<double,3> scalingX_;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
};

}

}

#endif // GINLA_EPETRAFVM_EVALUATOR_H
