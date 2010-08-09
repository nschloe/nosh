/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Sch\"omer

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

#include "Ginla_FVM_ModelEvaluator.h"

#include "Ginla_Komplex_LinearProblem.h"
#include "Ginla_Komplex_DoubleMatrix.h"
#include "Ginla_FVM_State.h"
#include "Ginla_MagneticVectorPotential_Centered.h"

#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

// ============================================================================
Ginla::FVM::ModelEvaluator::
ModelEvaluator ( const Teuchos::RCP<VIO::Mesh::Mesh>                          & mesh,
                 const Teuchos::ParameterList                                 & params,
                 const Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> & mvp,
                 const Teuchos::RCP<Ginla::Komplex::LinearProblem>            & komplex
               ) :
        komplex_ ( komplex ),
        mvp_( mvp ),
        x_(  Teuchos::null ),
        firstTime_ ( true ),
        numParams_( 1 ),
        p_map_( Teuchos::null ),
        p_init_( Teuchos::null ),
        p_names_( Teuchos::null ),
        mesh_( mesh ),
        kineticEnergyOperatorGraph_( Teuchos::null ),
        kineticEnergyOperator_( Teuchos::null ),
        dKineticEnergyDMuOperator_( Teuchos::null ),
        kineticEnergyOperatorsMu_( 0.0 ),
        jacobianOperator_( Teuchos::rcp( new Ginla::Komplex::DoubleMatrix( komplex_->getComplexMap(),
                                                                           komplex_->getComplexMap()
                                                                         )
                                       )
                         ),
        graph_( Teuchos::null )
{ 
  mesh_->computeFvmEntities();
  
  this->setupParameters_( params );
  
  // compute stiffness matrix and load vector
  this->createKineticEnergyOperatorGraph_();

  kineticEnergyOperator_     = Teuchos::rcp( new ComplexMatrix( kineticEnergyOperatorGraph_ ) );
  dKineticEnergyDMuOperator_ = Teuchos::rcp( new ComplexMatrix( kineticEnergyOperatorGraph_ ) );
  
  this->assembleKineticEnergyOperators_( 0.0 );
  
  // prepare initial guess
  Teuchos::RCP<Ginla::FVM::State> initialState =
      Teuchos::rcp( new Ginla::FVM::State( komplex_->getComplexMap(), mesh_) );
  initialState->getPsiNonConst()->putScalar( double_complex(1.0,0.0) );
//   initialState->getPsiNonConst()->randomize();
  x_ = this->createSystemVector_( *initialState );
  
  return;
}
// ============================================================================
Ginla::FVM::ModelEvaluator::
~ModelEvaluator()
{
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
setupParameters_( const Teuchos::ParameterList & params )
{
  p_names_ = Teuchos::rcp( new Teuchos::Array<std::string>() );
  p_names_->append( "mu" );
  
  // setup parameter values
  numParams_ = p_names_->length();
  
  p_map_ = Teuchos::rcp( new Epetra_LocalMap( numParams_,
                                              0,
                                              komplex_->getRealMap()->Comm()
                                            )
                       );
  p_init_ = Teuchos::rcp(new Epetra_Vector(*p_map_));
  for ( int k=0; k<numParams_; k++ )
      (*p_init_)[k] = params.get<double>( (*p_names_)[k] );
  
  return;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FVM::ModelEvaluator::
get_x_map() const
{
  return komplex_->getRealMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FVM::ModelEvaluator::
get_f_map() const
{
  return komplex_->getRealMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::FVM::ModelEvaluator::
get_x_init () const
{
  TEUCHOS_ASSERT( !x_.is_null() );
  return x_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::FVM::ModelEvaluator::
get_p_init ( int l ) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_init_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FVM::ModelEvaluator::
get_p_map(int l) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_map_;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Ginla::FVM::ModelEvaluator::
get_p_names (int l) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_names_;
}
// ============================================================================
Teuchos::RCP<Epetra_Operator>
Ginla::FVM::ModelEvaluator::
create_W() const
{
  TEUCHOS_ASSERT( !komplex_.is_null() );
  return komplex_->getMatrix();
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
Ginla::FVM::ModelEvaluator::
createInArgs() const
{ 
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  
  inArgs.setModelEvalDescription( "FVM Ginzburg-Landau" );
  
  inArgs.set_Np( numParams_ );
  
  inArgs.setSupports( IN_ARG_x, true );
  
//   // for shifted matrix
//   inArgs.setSupports( IN_ARG_alpha, true );
//   inArgs.setSupports( IN_ARG_beta, true );
  
  return inArgs;
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
Ginla::FVM::ModelEvaluator::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  
  outArgs.setModelEvalDescription( "FVM Ginzburg-Landau" );
  
  outArgs.set_Np_Ng( numParams_, 0 ); // return parameters p and solution g
  
  outArgs.setSupports( OUT_ARG_f, true );
  outArgs.setSupports( OUT_ARG_DfDp, 0, DerivativeSupport(DERIV_MV_BY_COL) );
  outArgs.setSupports( OUT_ARG_W, true );
  
  outArgs.set_W_properties( DerivativeProperties( DERIV_LINEARITY_NONCONST,
                                                  DERIV_RANK_FULL,
                                                  false // supportsAdjoint
                                                )
                          );

  return outArgs;
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
evalModel( const InArgs  & inArgs, 
           const OutArgs & outArgs
         ) const
{
//   double alpha = inArgs.get_alpha();
//   double beta  = inArgs.get_beta();
  
  const Teuchos::RCP<const Epetra_Vector> & x_in = inArgs.get_x();
  

  // Parse InArgs
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  TEUCHOS_ASSERT( !p_in.is_null() );
  const double mu = (*p_in)[0];
  TEUCHOS_ASSERT( !isnan( mu ) );
  
  // compute F
  const Teuchos::RCP<Epetra_Vector>   f_out = outArgs.get_f();
  if ( !f_out.is_null() )
      this->computeF_( *x_in, mu, *f_out );

  // fill jacobian
  const Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  if( !W_out.is_null() )
  {
      Teuchos::RCP<Epetra_CrsMatrix> W_out_crs =
          Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>( W_out, true );
      this->computeJacobian_( *x_in, mu, *W_out_crs );
      
//       W_out_crs->Scale( beta );
//       
//       double diag = -alpha;
//       for ( int i=0; i<x_in->MyLength(); i++ )
//           W_out_crs->SumIntoMyValues( i, 1, &diag, &i );
  }

  // dF / dmu
  Teuchos::RCP<Epetra_MultiVector> dfdp_out;
  if ( outArgs.Np() > 0 )
  {
      dfdp_out = outArgs.get_DfDp(0).getMultiVector();
      if ( !dfdp_out.is_null() )
          this->computeDFDp_( *x_in, mu, *((*dfdp_out)(0)) );
  }

  return;
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
computeF_ ( const Epetra_Vector & x,
            const double          mu,
            Epetra_Vector       & FVec
          ) const
{ 
  // convert from x to state
  const Teuchos::RCP<Ginla::FVM::State> state = this->createState_( x );

  // compute the GL residual
  Teuchos::RCP<Ginla::FVM::State> res = 
      Teuchos::rcp( new Ginla::FVM::State( komplex_->getComplexMap()->getComm(),
                                      mesh_
                                    )
                  );
  
  // reassemble linear operator if necessary
  if ( kineticEnergyOperatorsMu_ != mu )
      this->assembleKineticEnergyOperators_( mu );

  // compute r = K*psi
  // TODO check alpha, beta
  kineticEnergyOperator_->multiply( *(state->getPsi()),
                                    *(res->getPsiNonConst()),
                                    Teuchos::NO_TRANS,
                                    double_complex( 1.0, 0.0 ),
                                    double_complex( 0.0, 0.0 )
                                  );

//   std::cout << "SSS" << std::endl;
//   Teuchos::ArrayRCP<const double_complex> rView = res->getPsi()->get1dView();
//   for ( int k=0; k<rView.size(); k++ )
//     std::cout << k << " " << rView[k] << std::endl;

  
  // add nonlinear part (mass lumping)
  Teuchos::ArrayRCP<const double_complex> psiView = state->getPsi()->get1dView();
  Teuchos::ArrayRCP<double_complex>       resView = res->getPsiNonConst()->get1dViewNonConst();
  Teuchos::ArrayRCP<const double>         cvView  = mesh_->getControlVolumes()->get1dView();
  for ( int k=0; k<cvView.size(); k++ )
      resView[k] -= cvView[k] * psiView[k] * ( 1.0 - std::norm(psiView[k]) );
  
  // TODO Avoid this explicit copy?
  // transform back to fully real equation
  FVec = *(this->createSystemVector_( *res ));

  

//   TEUCHOS_ASSERT( !fvmLaplacian_.is_null() );
//   
//   FVec.PutScalar( 0.0 );
//   
//   // convert the element to a complex-valued vector
//   psi = ;
//   
//   // compute r = M*x
//   fvmLaplacian_->Multiply( false, x, FVec );
// 
//   // add nonlinear part to r.
//   // mass lumping  
//   Teuchos::ArrayRCP<Point> nodes = mesh_->getNodes();
//   Teuchos::RCP<DoubleVector> controlVolumes = mesh_->getControlVolumes();
//   Teuchos::ArrayRCP<const double> cvView = controlVolumes->get1dView();
//   for ( int k=0; k<nodes.size(); k++ )
//   {
//       double val =  mu * exp( x[k] ) * cvView[k]; // mass lumping
//       FVec.SumIntoGlobalValues( 1, &val, &k );
//   }

  return;
}
// ============================================================================
Teuchos::RCP<Ginla::FVM::State>
Ginla::FVM::ModelEvaluator::
createState_( const Epetra_Vector & x ) const
{
    const Teuchos::RCP<ComplexVector> psi = komplex_->real2complex ( x );
    return Teuchos::rcp( new Ginla::FVM::State( psi, mesh_ ) );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::FVM::ModelEvaluator::
createSystemVector_(  const Ginla::FVM::State & state ) const
{
    return komplex_->complex2real ( state.getPsi() );
}
// =============================================================================
void
Ginla::FVM::ModelEvaluator::
assembleKineticEnergyOperators_( const double mu ) const
{ 
  kineticEnergyOperator_->setAllToScalar( 0.0 );
  dKineticEnergyDMuOperator_->setAllToScalar( 0.0 );

  mvp_->setH0( mu );
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the elements, create local load vector and mass matrix,
  // and insert them into the global matrix.
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems = mesh_->getElems();
  Teuchos::ArrayRCP<Point> nodes = mesh_->getNodes();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeLengths = mesh_->getEdgeLengths();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > coedgeLengths = mesh_->getCoedgeLengths();
  
  TEUCHOS_ASSERT( !elems.is_null() );
  TEUCHOS_ASSERT( !nodes.is_null() );
  TEUCHOS_ASSERT( !edgeLengths.is_null() );
  TEUCHOS_ASSERT( !coedgeLengths.is_null() );

  for ( int k=0; k<elems.size(); k++ )
  {
      Teuchos::ArrayRCP<ORD> & elem = elems[k];
      // loop over the edges
      int n = elem.size();
      for ( int l=0; l<n; l++ )
      { 
          // two subsequent indices determine an edge
          Teuchos::Tuple<ORD,2> nodeIndices;
          nodeIndices[0] = elem[l];
          nodeIndices[1] = elem[(l+1)%n];
          
          // edge / co-edge ratio
          double alpha = coedgeLengths[k][l]
                       / edgeLengths  [k][l];

          // -------------------------------------------------------------------
          // Compute the integral
          //
          //    I = \int_{x0}^{xj} (xj-x0).A(x) dx
          //
          // numerically by the midpoint rule, i.e.,
          //
          //    I ~ |xj-x0| * (xj-x0) . A( 0.5*(xj+x0) ).
          //
          Point midpoint; // get A(midpoint)
          Point & node0 = nodes[ nodeIndices[0] ];
          Point & node1 = nodes[ nodeIndices[1] ];
          for (int i=0; i<midpoint.size(); i++ )
              midpoint[i] = 0.5 * ( node0[i] + node1[i] );

          Teuchos::Tuple<double_complex,2> values;
          // -------------------------------------------------------------------
          Teuchos::RCP<Point> a = mvp_->getA( midpoint );

          // project it onto the edge
          double aProj = 0.0;
          for (int i=0; i<midpoint.size(); i++ )
              aProj += ( node1[i] - node0[i] ) * (*a)[i];
          
          const double aInt = aProj * edgeLengths[k][l];

          values[0] =   alpha;
          values[1] = - alpha * exp( -IM * aInt );
          kineticEnergyOperator_->sumIntoGlobalValues( nodeIndices[0], nodeIndices, values );
                                  
          values[0] = - alpha * exp(  IM * aInt ); // integration the other way around
          values[1] =   alpha;
          kineticEnergyOperator_->sumIntoGlobalValues( nodeIndices[1], nodeIndices, values );
          
          // -------------------------------------------------------------------
          // Derivative by \mu.
          // Project derivative onto the edge.
          Teuchos::RCP<Point> Da = mvp_->getDADh0( midpoint );
          double DaProj = 0.0;
          for (int i=0; i<midpoint.size(); i++ )
              DaProj += ( node1[i] - node0[i] ) * (*Da)[i];
          
          const double DaInt = DaProj * edgeLengths[k][l];

          values[0] =   alpha;
          values[1] = - alpha * (-IM * DaInt ) * exp( -IM * aInt );
          dKineticEnergyDMuOperator_->sumIntoGlobalValues( nodeIndices[0], nodeIndices, values );
          
          values[0] = - alpha * ( IM * DaInt ) * exp( IM * aInt );
          values[1] =   alpha;
          dKineticEnergyDMuOperator_->sumIntoGlobalValues( nodeIndices[0], nodeIndices, values );
          // -------------------------------------------------------------------
      }
  }  

  kineticEnergyOperator_->fillComplete();
  dKineticEnergyDMuOperator_->fillComplete();
  
  kineticEnergyOperatorsMu_ = mu;
  
  return;
}
// =============================================================================
void
Ginla::FVM::ModelEvaluator::
createKineticEnergyOperatorGraph_()
{
  TEUCHOS_ASSERT( !komplex_.is_null() );
    
  // allocate the graph
  // TODO don't allocate that huge a graph
  // This is a bug in Trilinos, wait unti this is fixed.
  int maxNumEntriesPerRow = komplex_->getComplexMap()->getGlobalNumElements();
  kineticEnergyOperatorGraph_ = Teuchos::rcp( new TCrsGraph( komplex_->getComplexMap(),
                                                             komplex_->getComplexMap(),
                                                             maxNumEntriesPerRow
                                                           )
                                            );

  TEUCHOS_ASSERT( !mesh_.is_null() );
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems = mesh_->getElems();
  
  // connect all nodes with all other nodes in each element
  int n = elems.size();
  for ( int k=0; k<n; k++ )
  {
      Teuchos::ArrayRCP<ORD> & elem = elems[k];
      for ( int l=0; l<elem.size(); l++ )
          kineticEnergyOperatorGraph_->insertLocalIndices( elem[l], elem() );
  }

  kineticEnergyOperatorGraph_->fillComplete();
  
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
computeDFDp_ ( const Epetra_Vector & x,
               const double          mu,
                     Epetra_Vector & FVec
             ) const
{
  FVec.PutScalar( 0.0 );

//   // reassemble linear operator if necessary
//   if ( kineticEnergyOperatorsMu_ != mu )
//       this->assembleKineticEnergyOperators_( mu );
//   
//   // convert from x to state
//   const Teuchos::RCP<const Ginla::FVM::State> state = this->createState_( x );
// 
//   // initialize the sensitivity
//   Teuchos::RCP<Ginla::FVM::State> sens = 
//       Teuchos::rcp( new Ginla::FVM::State( komplex_->getComplexMap()->getComm(),
//                                       mesh_
//                                     )
//                   );
// 
//   // compute r = K*psi
//   dKineticEnergyDMuOperator_->multiply( *(state->getPsi()),
//                                         *(sens->getPsiNonConst()),
//                                         Teuchos::NO_TRANS,
//                                         double_complex( 1.0, 0.0 ),
//                                         double_complex( 0.0, 0.0 )
//                                       );
//   
//   // TODO Avoid this explicit copy?
//   // transform back to fully real equation
//   FVec = *(this->createSystemVector_( *sens ));

  return;
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
computeJacobian_ ( const Epetra_Vector & x,
                   const double          mu,
                   Epetra_CrsMatrix    & Jac
                 ) const
{
  // reassemble linear operator if necessary
  if ( kineticEnergyOperatorsMu_ != mu )
      this->assembleKineticEnergyOperators_( mu );
  
  // convert from x to state
  const Teuchos::RCP<const Ginla::FVM::State> state = this->createState_( x );

  // construct parts A and B of the linear operator part
  
  // A = K + I * ( 1 - 2*|psi|^2 )
  // B = diag( -psi^2 )
//   Teuchos::RCP<ComplexMatrix> A = this->deepCopy_( kineticEnergyOperator_ );
//   Teuchos::RCP<ComplexMatrix> B =
//       Teuchos::rcp( new ComplexMatrix( komplex_->getComplexMap(), 1, Tpetra::StaticProfile ) );
  Teuchos::ArrayRCP<const double_complex> psiView = state->getPsi()->get1dView();
  Teuchos::ArrayRCP<const double>         cvView  = mesh_->getControlVolumes()->get1dView();
  // update the diagonal values
  Teuchos::ArrayRCP<const ORD> indices;
  Teuchos::ArrayRCP<const double_complex> values;
  for ( unsigned int k=0; k<psiView.size(); k++ )
  {
      // get the local row view of the kinetic energy operator
      kineticEnergyOperator_->getLocalRowView( k, indices, values );
      jacobianOperator_->putALocalValues( k, indices, values );
      
      // Update the diagonal.
      // Unfortunately we have to artifically convert alpha back to
      // complex-valued as otherweise sumIntoGlobalValues won't eat it.
      double_complex alpha = - (1.0 - 2.0*norm(psiView[k])) * cvView[k];
      jacobianOperator_->AsumIntoGlobalValues( k, Teuchos::tuple<ORD>(k), Teuchos::tuple(alpha) );
      
      double_complex beta = psiView[k]*psiView[k] * cvView[k];
      jacobianOperator_->putBLocalValues( k, Teuchos::tuple<ORD>(k), Teuchos::tuple( beta ) );
  }

  // finalize the Jacobian operator
  if ( firstTime_ )
    jacobianOperator_->finalize();

  // feed it into komplex_
  komplex_->zeroOutMatrix();
  komplex_->update( jacobianOperator_, firstTime_ );
  
  if ( firstTime_ )
  {
      komplex_->finalizeMatrix();
      firstTime_ = false;
  }
  
  // TODO
  // avoid this explicit copy?
  Jac = *(komplex_->getMatrix());

  return;
}
// =============================================================================
Teuchos::RCP<ComplexMatrix>
Ginla::FVM::ModelEvaluator::
deepCopy_ ( const Teuchos::RCP<const ComplexMatrix> & A
          ) const
{
  Teuchos::RCP<ComplexMatrix> B =
      Teuchos::rcp( new ComplexMatrix( A->getCrsGraph() ) );
      
  // copy row by row
  Teuchos::ArrayRCP<const ORD> indices;
  Teuchos::ArrayRCP<const double_complex> values;
  for ( ORD k; k<A->getRowMap()->getNodeNumElements(); k++ )
  {
      A->getLocalRowView( k, indices, values );
      B->replaceLocalValues( k, indices(), values() );
  }
  
  return B;
}
// =============================================================================