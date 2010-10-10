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

#ifndef GINLA_FVM_MODELEVALUATOR_H
#define GINLA_FVM_MODELEVALUATOR_H
// -----------------------------------------------------------------------------
// includes
#include <EpetraExt_ModelEvaluator.h>

#include "VIO_Mesh_Mesh.h"
#include "Ginla_StateTranslator.h"
#include "Ginla_ParameterHost_Virtual.h"

#include <Teuchos_Comm.hpp>
#include <Teuchos_map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
// -----------------------------------------------------------------------------
typedef Tpetra::CrsGraph<ORD> TCrsGraph;
typedef Tpetra::CrsMatrix<double_complex,ORD> ComplexMatrix;
// -----------------------------------------------------------------------------
// forward declarations
namespace Ginla
{
  namespace MagneticVectorPotential
  {
    class Virtual;
  }
  namespace FVM {
    class State;
  }
}
namespace Komplex2
{
  class LinearProblem;
  class DoubleMatrix;
}

class Epetra_CrsGraph;
// -----------------------------------------------------------------------------
namespace Ginla {
namespace FVM {
  
class ModelEvaluator:
    public EpetraExt::ModelEvaluator,
    public Ginla::StateTranslator,
    public Ginla::ParameterHost::Virtual
{

public:
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  //! Constructor without initial guess.
  ModelEvaluator ( const Teuchos::RCP<VIO::Mesh::Mesh>                         & mesh,
                   const Teuchos::ParameterList                                & params,
                   const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp,
                   const Teuchos::RCP<Komplex2::LinearProblem>                 & komplex,
                   const Teuchos::RCP<Ginla::FVM::State>                       & initialState
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
    Teuchos::RCP<Ginla::State::Virtual>
    createState( const Epetra_Vector & x ) const;
    
    Teuchos::RCP<Epetra_Vector>
    createSystemVector(  const Ginla::State::Virtual & state ) const;

    void
    createSystemVector( const Ginla::State::Virtual & state,
                              Epetra_Vector         & x
                      ) const;
    
    virtual
    Teuchos::RCP<LOCA::ParameterVector>
    getParameters() const;
  
private:
  void
  computeF_ ( const Epetra_Vector            & x,
              const double                     lambda,
              const Teuchos::Tuple<double,3> & scaling,
                    Epetra_Vector            & FVec
            ) const;
              
  void
  computeDFDp_ ( const Epetra_Vector            & x,
                 const double                     mu,
                 const Teuchos::Tuple<double,3> & scaling,
                       Epetra_Vector            & FVec
               ) const;
              
  void
  computeJacobian_ ( const Epetra_Vector            & x,
                     const double                     lambda,
                     const Teuchos::Tuple<double,3> & scaling,
                           Epetra_CrsMatrix         & Jac
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

   const Teuchos::RCP<Komplex2::LinearProblem> komplex_;
   const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp_;

   Teuchos::RCP<Epetra_Vector> x_;
   mutable bool firstTime_;
   
   int numParams_;
   
   Teuchos::RCP<Epetra_Map> p_map_;
   Teuchos::RCP<Epetra_Vector> p_init_;
   Teuchos::RCP<Teuchos::Array<std::string> > p_names_;
   
   const Teuchos::RCP<VIO::Mesh::Mesh> mesh_;
   Teuchos::RCP<TCrsGraph> kineticEnergyOperatorGraph_;
   mutable Teuchos::RCP<ComplexMatrix> kineticEnergyOperator_;
   mutable Teuchos::RCP<ComplexMatrix> dKineticEnergyDMuOperator_;
   mutable bool kineticEnergyOperatorsAssembled_;
   mutable double kineticEnergyOperatorsMu_;
   mutable Teuchos::Tuple<double,3> kineticEnergyOperatorsScaling_;
   Teuchos::RCP<Komplex2::DoubleMatrix> jacobianOperator_;
//    Teuchos::RCP<Epetra_CrsMatrix> fvmLaplacian_;
//    Teuchos::RCP<Epetra_CrsMatrix> jacobian_;
   Teuchos::RCP<Epetra_CrsGraph> graph_;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
private:  
    void
    setupParameters_( const Teuchos::ParameterList & params );
    
    void
    assembleKineticEnergyOperators_( double                           mu,
                                     const Teuchos::Tuple<double,3> & scaling
                                   ) const;

    void
    createKineticEnergyOperatorGraph_();
    
    Teuchos::RCP<ComplexMatrix>
    deepCopy_ ( const Teuchos::RCP<const ComplexMatrix> & A
              ) const;
              
    bool
    kineticEnergyOperatorsUpToDate_( const double                     mu,
                                     const Teuchos::Tuple<double,3> & scaling
                                   ) const;
              
    // this value is really just for getParameters()
    mutable double mu_;
    mutable double scaling_;
    mutable Teuchos::Tuple<double,3> scalingX_;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
};

}

}

#endif // GINLA_FVM_MODELEVALUATOR_H
