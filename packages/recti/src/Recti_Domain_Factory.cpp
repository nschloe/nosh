/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Sch\"omer

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

#include "Recti_Domain_Factory.h"

#include "Recti_Domain_Circle.h"
#include "Recti_Domain_Ellipse.h"
#include "Recti_Domain_Polygon.h"
#include "Recti_Domain_Rectangle.h"
#include "Recti_Domain_Square.h"

// ============================================================================
Recti::Domain::Factory::
Factory( const Teuchos::ParameterList & pList ):
    pList_( pList )
{ 
}
// ============================================================================    
Recti::Domain::Factory::
~Factory()
{
}
// ============================================================================
Teuchos::RCP<Recti::Domain::Abstract>
Recti::Domain::Factory::
build() const
{
  std::string domainType = pList_.get<std::string>("Type");
  if ( domainType.compare("Circle")==0 )
      return buildCircle();
  else if ( domainType.compare("Ellipse")==0 )
      return buildEllipse();
  else if ( domainType.compare("Polygon")==0 )
      return buildPolygon();
  else if ( domainType.compare("Rectangle")==0 )
      return buildRectangle();
  else if ( domainType.compare("Square")==0 )
      return buildSquare();
  else
  {
      TEST_FOR_EXCEPTION( true, std::logic_error,
                          "Unknown domain type \"" << domainType << "\"." );
  }
}
// ============================================================================
Teuchos::RCP<Recti::Domain::Abstract>
Recti::Domain::Factory::
buildCircle() const
{
    double radius = pList_.get<double>("Radius"); 
    return Teuchos::rcp( new Circle( radius ) );
}
// ============================================================================
Teuchos::RCP<Recti::Domain::Abstract>
Recti::Domain::Factory::
buildEllipse() const
{
    double tDiam = pList_.get<double>("Transverse diameter");
    double cDiam = pList_.get<double>("Conjugate diameter");
    return Teuchos::rcp( new Ellipse( tDiam, cDiam ) );
}
// ============================================================================
Teuchos::RCP<Recti::Domain::Abstract>
Recti::Domain::Factory::
buildPolygon() const
{
    const Teuchos::ParameterList & pointsList =
            pList_.sublist("Vertices");
    
    Teuchos::Array<DoubleTuple> vertices;
    
    Teuchos::ParameterList::ConstIterator k;
    for ( k=pointsList.begin(); k!=pointsList.end(); ++k )
    {
        Teuchos::ParameterEntry e = pointsList.entry(k);
        TEUCHOS_ASSERT( e.isList() );
        
        std::string label = pointsList.name(k);
        const Teuchos::ParameterList & coordinates =
                pointsList.sublist(label);
        
        DoubleTuple X = Teuchos::tuple( coordinates.get<double>("x"),
                                        coordinates.get<double>("y") );
        vertices.push_back( X );
    }
    
    return Teuchos::rcp( new Polygon( vertices ) );
}
// ============================================================================
Teuchos::RCP<Recti::Domain::Abstract>
Recti::Domain::Factory::
buildRectangle() const
{
    double width = pList_.get<double>("Width");
    double height = pList_.get<double>("Height");

    return Teuchos::rcp( new Rectangle( width, height ) );
}
// ============================================================================
Teuchos::RCP<Recti::Domain::Abstract>
Recti::Domain::Factory::
buildSquare() const
{
    double edgeLength = pList_.get<double>("Edge length"); 
    return Teuchos::rcp( new Square( edgeLength ) );
}
// ============================================================================
// ============================================================================
// nonmember functions
Teuchos::RCP<Recti::Domain::Abstract>
Recti::Domain::
buildDomain(  const Teuchos::ParameterList & pList )
{
  Recti::Domain::Factory factory( pList );
  return factory.build();
}
// ============================================================================