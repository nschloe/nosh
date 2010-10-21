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

#include "Ginla_FDM_ModelEvaluator_Bordered.h"

#include <Epetra_LocalMap.h>
#include <Epetra_CrsMatrix.h>

#include "Ginla_FDM_Operator_Virtual.h"
#include "Komplex2_LinearProblem.h"

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Comm.h>

#include <AztecOO_ConditionNumber.h>

// ============================================================================
Ginla::FDM::ModelEvaluator::Bordered::
Bordered ( const Teuchos::RCP<Ginla::FDM::Operator::Virtual> & glOperator,
           const Teuchos::RCP<Komplex2::LinearProblem>       & komplex,
           const Teuchos::RCP<Recti::Grid::General>          & grid,
           const Teuchos::ParameterList                      & params
         ) :
         modelEvaluatorDefault_ ( Teuchos::rcp( new Ginla::FDM::ModelEvaluator::Default( glOperator, komplex, grid, params ) ) ),
         firstTime_ ( true ),
         regularMap_ (  modelEvaluatorDefault_->getMap() ),
         extendedMap_ ( this->createExtendedMap_ ( *regularMap_ ) ),
         importFromExtendedMap_(  Epetra_Import(*regularMap_,*extendedMap_) ),
         importFromRegularMap_( Epetra_Import(*extendedMap_,*regularMap_) ),
         jacobian_ ( new Epetra_CrsMatrix ( Copy, *extendedMap_, 0 ) ),
         x_(  Teuchos::null )
{ 
  Teuchos::RCP<Ginla::State::Updatable> initialState =
      Teuchos::rcp( new Ginla::FDM::State( komplex->getComplexMap()->getComm(),
                                           grid ) );
                                      
  TEUCHOS_ASSERT( initialState->getPsi()->getMap()->isSameAs( *komplex->getComplexMap()) );
                                      
  initialState->getPsiNonConst()->putScalar( double_complex(1.0,0.0) );
  
  *x_ = *(this->createSystemVector( *initialState ));
  
  return;
}
// ============================================================================
Ginla::FDM::ModelEvaluator::Bordered::
Bordered ( const Teuchos::RCP<Ginla::FDM::Operator::Virtual> & glOperator,
           const Teuchos::RCP<Komplex2::LinearProblem>       & komplex,
           const Ginla::State::Updatable                       & state,
           const Teuchos::ParameterList                      & params
         ) :
         modelEvaluatorDefault_ ( Teuchos::rcp( new Ginla::FDM::ModelEvaluator::Default( glOperator, komplex, state, params ) ) ),
         firstTime_ ( true ),
         regularMap_ (  modelEvaluatorDefault_->getMap() ),
         extendedMap_ ( this->createExtendedMap_ ( *regularMap_ ) ),
         importFromExtendedMap_(  Epetra_Import(*regularMap_,*extendedMap_) ),
         importFromRegularMap_( Epetra_Import(*extendedMap_,*regularMap_) ),
         jacobian_ ( new Epetra_CrsMatrix ( Copy, *extendedMap_, 0 ) ),
         x_( this->createSystemVector( state ) )
{
  // make sure the maps are compatible
  TEUCHOS_ASSERT( state.getPsi()->getMap()->isSameAs( *komplex->getComplexMap()) );  
  return;
}
// ============================================================================
Ginla::FDM::ModelEvaluator::Bordered::
~Bordered()
{
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FDM::ModelEvaluator::Bordered::
get_x_map() const
{
  return extendedMap_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FDM::ModelEvaluator::Bordered::
get_f_map() const
{
  return extendedMap_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::FDM::ModelEvaluator::Bordered::
get_x_init () const
{
  return x_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::FDM::ModelEvaluator::Bordered::
get_p_init ( int l ) const
{
  return modelEvaluatorDefault_->get_p_init( l );
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FDM::ModelEvaluator::Bordered::
get_p_map(int l) const
{
  return modelEvaluatorDefault_->get_p_map( l );
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Ginla::FDM::ModelEvaluator::Bordered::
get_p_names (int l) const
{
  return modelEvaluatorDefault_->get_p_names( l );
}
// ============================================================================
Teuchos::RCP<Epetra_Operator>
Ginla::FDM::ModelEvaluator::Bordered::
create_W() const
{
  return jacobian_;
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
Ginla::FDM::ModelEvaluator::Bordered::
createInArgs() const
{
  return modelEvaluatorDefault_->createInArgs();
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
Ginla::FDM::ModelEvaluator::Bordered::
createOutArgs() const
{
  return modelEvaluatorDefault_->createOutArgs();
}
// ============================================================================
void
Ginla::FDM::ModelEvaluator::Bordered::
evalModel( const InArgs  & inArgs, 
           const OutArgs & outArgs
         ) const
{
  double alpha = inArgs.get_alpha();
  double beta  = inArgs.get_beta();
  
  const Teuchos::RCP<const Epetra_Vector> & x_in = inArgs.get_x();
  
  Teuchos::RCP<Epetra_Vector>   f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();

  Teuchos::RCP<Epetra_MultiVector> dfdp_out;
  if ( outArgs.Np() > 0 )
      dfdp_out = outArgs.get_DfDp(0).getMultiVector();
  
  // Parse InArgs
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  
  TEUCHOS_ASSERT( !p_in.is_null() );
  
  // build a LOCA::ParameterVector list to feed glOperator_.
  LOCA::ParameterVector locaParams;
  
  Teuchos::RCP<const Teuchos::Array<std::string> > p_names = modelEvaluatorDefault_->get_p_names( 0 );
  
  for ( int k=0; k<p_in->GlobalLength(); k++ )
      locaParams.addParameter( (*p_names)[k] , (*p_in)[k] ); 
  modelEvaluatorDefault_->getOperator()->setParameters( locaParams );
    
  // compute F
  if ( !f_out.is_null() )
      this->computeF_( *x_in, *f_out );

  // fill jacobian
  if( !W_out.is_null() )
  { 
      this->computeJacobian_( *x_in, *W_out );
      Teuchos::RCP<Epetra_CrsMatrix> W_out_crs =
          Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out, true);
      W_out_crs->Scale(beta);

      double diag = -alpha;
      for (int i=0; i<x_in->MyLength(); i++)
          W_out_crs->SumIntoMyValues(i, 1, &diag, &i);
  }
  
  // dF / dH_0
  if ( !dfdp_out.is_null() )
      this->computeDFDh0_( *x_in, *((*dfdp_out)(0)) );
  
  return;
}
// ============================================================================
void
Ginla::FDM::ModelEvaluator::Bordered::
computeF_ ( const Epetra_Vector & x,
            Epetra_Vector       & FVec
          ) const
{ 
    TEST_FOR_EXCEPTION ( regularMap_.is_null(),
                         std::logic_error,
                         "regularMap_ not properly initialized." );

    TEST_FOR_EXCEPTION ( extendedMap_.is_null(),
                         std::logic_error,
                         "extendedMap_ not properly initialized." );

    // make sure that the input and output vectors are correctly mapped
    TEST_FOR_EXCEPTION ( !x.Map().SameAs ( *extendedMap_ ),
                         std::logic_error,
                         "Maps of x and the computed real-valued map do not coincide. "
                         << "Check, for example, the number of elements "
                         << "(" << x.Map().NumGlobalElements() << " for x vs. "
                         << extendedMap_->NumGlobalElements() << " for extendedMap_)." );

    TEST_FOR_EXCEPTION ( !FVec.Map().SameAs ( *extendedMap_ ),
                         std::logic_error,
                         "Maps of FVec and the computed real-valued map do not coincide."
                         << "Check, for example, the number of elements "
                         << "(" << FVec.Map().NumGlobalElements() << " for FVec vs. "
                         << extendedMap_->NumGlobalElements() << " for extendedMap_)." );
    
    // strip off the phase constraint variable
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );

    double eta = x[ x.GlobalLength()-1 ];

    // compute GL in the default model
    Epetra_Vector shortFVec ( *regularMap_ );
    modelEvaluatorDefault_->computeF ( tmp, shortFVec );

    // \psi += -i\eta\psi
    Teuchos::RCP<Ginla::State::Updatable> state = modelEvaluatorDefault_->createState( tmp );
    state->getPsiNonConst()->scale( -IM*eta );
    TEUCHOS_ASSERT_EQUALITY( 0, tmp.Update( 1.0, *modelEvaluatorDefault_->createSystemVector( *state ), 1.0 ) );
    
    // copy over and add phase condition
    FVec.Import( shortFVec, importFromRegularMap_, Insert );

    // set last entry
    FVec.ReplaceGlobalValue ( shortFVec.GlobalLength(),
                              0,
                              0.0 );

    return;
}
// ============================================================================
void
Ginla::FDM::ModelEvaluator::Bordered::
computeDFDh0_ ( const Epetra_Vector & x,
                Epetra_Vector       & FVec
              ) const
{   
    TEST_FOR_EXCEPTION ( regularMap_.is_null(),
                         std::logic_error,
                         "regularMap_ not properly initialized." );

    TEST_FOR_EXCEPTION ( extendedMap_.is_null(),
                         std::logic_error,
                         "extendedMap_ not properly initialized." );

    // make sure that the input and output vectors are correctly mapped
    TEST_FOR_EXCEPTION ( !x.Map().SameAs ( *extendedMap_ ),
                         std::logic_error,
                         "Maps of x and the computed real-valued map do not coincide. "
                         << "Check, for example, the number of elements "
                         << "(" << x.Map().NumGlobalElements() << " for x vs. "
                         << extendedMap_->NumGlobalElements() << " for extendedMap_)." );

    TEST_FOR_EXCEPTION ( !FVec.Map().SameAs ( *extendedMap_ ),
                         std::logic_error,
                         "Maps of FVec and the computed real-valued map do not coincide."
                         << "Check, for example, the number of elements "
                         << "(" << FVec.Map().NumGlobalElements() << " for FVec vs. "
                         << extendedMap_->NumGlobalElements() << " for extendedMap_)." );
    
    // strip off the phase constraint variable
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );

    double eta = x[ x.GlobalLength()-1 ];

    // compute GL in the default model
    Epetra_Vector shortFVec ( *regularMap_ );
    modelEvaluatorDefault_->computeDFDh0 ( tmp, shortFVec );
    
    // copy over and add phase condition
    FVec.Import( shortFVec, importFromRegularMap_, Insert );

    // set last entry
    FVec.ReplaceGlobalValue ( shortFVec.GlobalLength(),
                              0,
                              0.0 );
    return;
}
// ============================================================================
void
Ginla::FDM::ModelEvaluator::Bordered::
computeJacobian_ ( const Epetra_Vector & x,
                   Epetra_Operator     & Jac
                 ) const
{ 
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );

    // Strip down Jac, too?
    // --   Not really necessary as it's not being used anyway.

    // compute the underlying Jacobian
    modelEvaluatorDefault_->computeJacobian ( tmp, Jac );
    
    TEST_FOR_EXCEPTION ( !extendedMap_.is_valid_ptr() || extendedMap_.is_null(),
                         std::logic_error,
                         "extendedMap_ not properly initialized." );
    
    double eta = x[ x.GlobalLength()-1 ];

    // TODO don't explicitly construct psi? get1dCopy on the rhs
    Teuchos::RCP<Ginla::State::Updatable> state = modelEvaluatorDefault_->createState( tmp );
    
    Teuchos::RCP<const ComplexVector> psi = state->getPsi();

    // get the unbordered Jacobian
    Teuchos::RCP<Epetra_CrsMatrix> regularJacobian = modelEvaluatorDefault_->getJacobianNonConst();
    
    // add "-i*eta" on the diagonal
    Epetra_Vector newDiag ( tmp );
    regularJacobian->ExtractDiagonalCopy( newDiag );
    Teuchos::RCP<Ginla::State::Updatable> diag = modelEvaluatorDefault_->createState( newDiag );
    diag->getPsiNonConst()->putScalar( -IM*eta );
    newDiag.Update( 1.0, *modelEvaluatorDefault_->createSystemVector(*diag), 1.0 );
    regularJacobian->ReplaceDiagonalValues( newDiag );
    
    // TODO: Conversion to real-valued vector in one go?
    // right bordering: (phi:=) -i*psi
    Teuchos::RCP<Ginla::FDM::State> stateFdm =
        Teuchos::rcp_dynamic_cast<Ginla::FDM::State>( state );
    Ginla::FDM::State phiState = *stateFdm;
    phiState.getPsiNonConst()->scale ( -IM );
    Teuchos::RCP<Epetra_Vector> rightBorder = modelEvaluatorDefault_->createSystemVector( phiState );

    // Get the lower bordering
    //     
    //    Im( psi_{old}^H, dpsi ) 
    //  = \Im( \sum_{k=1}^n \psi^*_k dpsi_k )
    //  = \sum_{k=1}^n ( -Im(\psi_k) Re(dpsi_k) + Re(\psi_k) Im(dpsi_k)
    //
    // Hence, the bordering tranlates into
    //
    //  lb = IM*psi_{old}
    //
    // in the complex state notation.
    Ginla::FDM::State lbState = *stateFdm;
    lbState.getPsiNonConst()->scale( IM );
    Teuchos::RCP<Epetra_Vector> lowerBorder = modelEvaluatorDefault_->createSystemVector( lbState );
    
    // corner element
    double d = 0.0;

    // create the bordered Jacobian out of this
    this->fillBorderedMatrix_ ( jacobian_,
                                regularJacobian,
                                *rightBorder,
                                *lowerBorder,
                                d,
                                firstTime_ );

    if ( firstTime_ )
    {
        TEUCHOS_ASSERT_EQUALITY ( 0, jacobian_->FillComplete() );
        TEUCHOS_ASSERT_EQUALITY ( 0, jacobian_->OptimizeStorage() );
        firstTime_ = false;
    }

    return;
}
// ============================================================================
Teuchos::RCP<Ginla::State::Updatable>
Ginla::FDM::ModelEvaluator::Bordered::
createState( const Epetra_Vector & x ) const
{
  // TODO Use chi as well!
  
  Epetra_Vector tmp ( *regularMap_ );
  tmp.Import( x, importFromExtendedMap_, Insert );
  
  return modelEvaluatorDefault_->createState( tmp );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::FDM::ModelEvaluator::Bordered::
createSystemVector( const Ginla::State::Updatable & state ) const
{
    TEUCHOS_ASSERT( !extendedMap_.is_null() );
    Teuchos::RCP<Epetra_Vector> x =
        Teuchos::rcp( new Epetra_Vector( *extendedMap_ ) );

    this->createSystemVector( state, *x );
                
    return x;
}
// =============================================================================
void
Ginla::FDM::ModelEvaluator::Bordered::
createSystemVector( const Ginla::State::Updatable & state,
                          Epetra_Vector         & x
                  ) const
{
    TEUCHOS_ASSERT( x.Map().SameAs( *extendedMap_ ) );
    TEUCHOS_ASSERT( !modelEvaluatorDefault_.is_null() );

    x.Import( *modelEvaluatorDefault_->createSystemVector( state ),
              importFromRegularMap_,
              Insert );

    // set last entry
    x.ReplaceGlobalValue ( x.GlobalLength()-1, 0, state.getChi() );
                
    return;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Ginla::FDM::ModelEvaluator::Bordered::
createExtendedMap_ ( const Epetra_BlockMap & realMap ) const
{
    // fill up realMapGIDs
    int numMyElements = realMap.NumMyElements();
    Teuchos::Array<int> myElements ( numMyElements );
    realMap.MyGlobalElements ( myElements.getRawPtr() );

    // Unconditionally put the phase constraint on the first process.
    int myPID = realMap.Comm().MyPID();
    if ( myPID==0 )
    {
        int n = realMap.NumGlobalElements();
        // extend the GIDs by the phase constraint
        myElements.append ( n );
    }

    int numGlobalElements = realMap.NumGlobalElements() + 1;
    return Teuchos::rcp ( new Epetra_Map ( numGlobalElements,
                                           myElements.length(),
                                           myElements.getRawPtr(),
                                           realMap.IndexBase(),
                                           realMap.Comm() )
                        );
}
// =============================================================================
void
Ginla::FDM::ModelEvaluator::Bordered::
fillBorderedMatrix_ ( const Teuchos::RCP<      Epetra_CrsMatrix> & extendedMatrix,
                      const Teuchos::RCP<const Epetra_CrsMatrix> & regularMatrix,
                      const Epetra_Vector                        & rightBorder,
                      // TODO Declare the following const as soon as Trilinos allows (ReplaceGlobalValues)
                      Epetra_Vector                              & lowerBorder,
                      double                                       d,
                      bool                                         firstTime
                    ) const
{
    TEUCHOS_ASSERT ( regularMatrix.is_valid_ptr()  && !regularMatrix.is_null() );
    TEUCHOS_ASSERT ( extendedMatrix.is_valid_ptr() && !extendedMatrix.is_null() );

    int m = regularMatrix->NumGlobalRows();
    int n = regularMatrix->NumGlobalCols();

    // check if the sizes all match
    TEUCHOS_ASSERT_EQUALITY ( m+1, extendedMatrix->NumGlobalRows() );
    TEUCHOS_ASSERT_EQUALITY ( n+1, extendedMatrix->NumGlobalCols() );

    // make sure the maps coincide
    TEUCHOS_ASSERT ( lowerBorder.Map().SameAs ( regularMatrix->OperatorDomainMap() ) );
    TEUCHOS_ASSERT ( rightBorder.Map().SameAs ( regularMatrix->OperatorRangeMap() ) );

    int numMyRows = regularMatrix->NumMyRows();

    // fill the matrix with the entries
    int numRowNonZeros;

    int maxNumEntries = regularMatrix->MaxNumEntries() + 1; // count the last column in
    int    * indices = new int   [maxNumEntries];
    double * values  = new double[maxNumEntries];
    for ( int myRow=0; myRow<numMyRows; myRow++ )
    {
        // extract row view
        TEUCHOS_ASSERT_EQUALITY ( 0, regularMatrix->ExtractMyRowView ( myRow, numRowNonZeros, values, indices ) );

        // *indices are given in global indexing.
        int globalRow = extendedMatrix->Map().GID ( myRow );

        // Write the data to the new matrix.
        // Only panic for negative return codes.
        TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutRow_ ( extendedMatrix, globalRow, numRowNonZeros, values, indices, firstTime ) );

        // add last column
        double val = rightBorder[globalRow];
        TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutRow_ ( extendedMatrix, globalRow, 1, &val, &n, firstTime ) );
    }


    // set the last row
    // put the last row piece by piece
    int    numMyElements = lowerBorder.Map().NumMyElements();
    std::vector<double> myLowerBorderValues ( numMyElements );
    // TODO only use ExtractView
    lowerBorder.ExtractCopy ( & ( myLowerBorderValues[0] ) );
    int * myLowerBorderIndices = new int[numMyElements];
    for ( int k=0; k<numMyElements; k++ )
        myLowerBorderIndices[k] = lowerBorder.Map().GID ( k );
    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutRow_ ( extendedMatrix, n, numMyElements, & ( myLowerBorderValues[0] ), myLowerBorderIndices, firstTime ) );

    // set the last element d
    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutRow_ ( extendedMatrix, n, 1, &d, &n, firstTime ) );

    return;
}
// =============================================================================
int
Ginla::FDM::ModelEvaluator::Bordered::
PutRow_ ( const Teuchos::RCP<Epetra_CrsMatrix> A,
          const int                            globalRow,
          const int                            numIndices,
          double                             * values,
          int                                * indices,
          const bool                           firstTime
        ) const
{
    if ( firstTime )
        return A->InsertGlobalValues ( globalRow, numIndices, values, indices );
    else
        return A->ReplaceGlobalValues ( globalRow, numIndices, values, indices );
}
// =============================================================================
