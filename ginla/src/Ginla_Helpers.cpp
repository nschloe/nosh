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

#include "Ginla_Helpers.hpp"

#include "NOX_Abstract_Vector.H"
#include "NOX_Epetra_Vector.H"

#include <Tpetra_Vector.hpp>

typedef std::complex<double> double_complex;

namespace Ginla {
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Helpers::
teuchosParameterList2locaParameterVector( const Teuchos::ParameterList & p
                                        )
{
  Teuchos::RCP<LOCA::ParameterVector> pL =
      Teuchos::rcp( new LOCA::ParameterVector() );
        
  Teuchos::ParameterList::ConstIterator k;
  double * dummy = NULL;
  for ( k=p.begin(); k!=p.end(); ++k )
  {
    Teuchos::ParameterEntry e = p.entry(k);
    if ( e.isType<double>() )
      pL->addParameter( p.name(k),
                        e.getValue<double>(dummy) );
  }
  
  return pL;
}
// ============================================================================
Teuchos::RCP<Teuchos::ParameterList>
Helpers::
locaParameterVector2teuchosParameterList( const LOCA::ParameterVector & pL )
{
    Teuchos::RCP<Teuchos::ParameterList> p =
      Teuchos::rcp( new Teuchos::ParameterList() );
      
    appendToTeuchosParameterList(  *p, pL );

    return p;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Helpers::
mergeLocaParameterVectors( const LOCA::ParameterVector & p0,
                           const LOCA::ParameterVector & p1
                         )
{
  // intialize p with p0
  Teuchos::RCP<LOCA::ParameterVector> p = 
          Teuchos::rcp( new LOCA::ParameterVector( p0 ) );

  // add elements from p1
  for( int k=0; k<p1.length(); k++ )
  {
    double value = p1.getValue(k);
    std::string label = p1.getLabel(k);
    if( p->isParameter(label) )
    {
        // If the entry already exists, make sure the values
        // coincide.
        TEUCHOS_ASSERT_EQUALITY( p->getValue(label), value );
    }
    else
    {
        p->addParameter( label, value );
    }
  }
  
  return p;
}
// ============================================================================
void
Helpers::
appendToTeuchosParameterList( Teuchos::ParameterList      & p,
                              const LOCA::ParameterVector & pL,
                              const std::string           & labelPrepend
                            )
{       
  Teuchos::ParameterList::ConstIterator k;
  for ( int k=0; k<pL.length(); k++ )
     p.set<double>( labelPrepend + pL.getLabel(k), pL[k] );  
  
  return;
}
// ============================================================================
unsigned int
Helpers::
numDigits ( const int i )
{
    int numDigits = 0;
    int ii = i;
    if ( ii < 0 )
        ii = -ii;

    while ( ii > 0 )
    {
        numDigits++;
        ii/=10;
    }
    return numDigits;
}
// ============================================================================
} // namespace Ginla