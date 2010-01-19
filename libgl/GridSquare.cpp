/*
 * Grid.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schl�mer
 */

#include "GridSquare.h"

#include "ioVirtual.h"
#include "ioFactory.h"

// =============================================================================
// Class constructor
GridSquare::GridSquare(Teuchos::Tuple<unsigned int,2> Nx, double scaling) :
  GridVirtual( scaling,
               Teuchos::tuple<double>(scaling/Nx[0],scaling/Nx[1]),
               pow( scaling, 2 ),
               (Nx[0]+1)*(Nx[1]+1),
               2*( Nx[0]+Nx[1] )
             ),
  Nx_(Nx)
{
}
// =============================================================================
// Class constructor for this classes data only.
GridSquare::GridSquare( Teuchos::Tuple<unsigned int,2> Nx ) :
  GridVirtual(),
  Nx_(Nx)
{
}
// =============================================================================
// Destructor
GridSquare::~GridSquare()
{
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getX(Teuchos::Array<int> i) const
{
  Teuchos::RCP<Teuchos::Array<double> > x = Teuchos::rcp(new Teuchos::Array<
      double>(2));
  (*x)[0] = i[0] * h_[0];
  (*x)[1] = i[1] * h_[1];
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getXLeft(unsigned int k) const
{
  Teuchos::RCP<Teuchos::Array<int> >    i( k2i(k)   );
  Teuchos::RCP<Teuchos::Array<double> > x( getX(*i) );
  (*x)[0] -= 0.5 * h_[0];
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getXRight(unsigned int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[0] += 0.5 * h_[0];
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getXBelow(unsigned int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[1] -= 0.5 * h_[1];
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getXAbove(unsigned int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[1] += 0.5 * h_[1];
  return x;
}
// =============================================================================
unsigned int
GridSquare::getKLeft( unsigned int k ) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[0] -= 1;
  return i2k(i);
}
// =============================================================================
unsigned int
GridSquare::getKRight(unsigned int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[0] += 1;
  return i2k(i);
}
// =============================================================================
unsigned int
GridSquare::getKBelow(unsigned int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[1] -= 1;
  return i2k(i);
}
// =============================================================================
unsigned int
GridSquare::getKAbove(unsigned int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[1] += 1;
  return i2k(i);
}
// =============================================================================
// maps a running index k to a 2D index i
Teuchos::RCP<Teuchos::Array<int> >
GridSquare::k2i(unsigned int k) const
{
  int d = 2;
  Teuchos::RCP<Teuchos::Array<int> > i = Teuchos::rcp(
      new Teuchos::Array<int>(d));

  if (k < Nx_[0])
    { // south
      (*i)[0] = k;
      (*i)[1] = 0;
    }
  else if (k < Nx_[0]+Nx_[1])
    { // east
      (*i)[0] = Nx_[0];
      (*i)[1] = k - Nx_[0];
    }
  else if (k < 2*Nx_[0]+Nx_[1])
    { // north
      (*i)[0] = 2*Nx_[0]+Nx_[1] -k;
      (*i)[1] = Nx_[1];
    }
  else if (k < 2*(Nx_[0]+Nx_[1]) )
    { // west
      (*i)[0] = 0;
      (*i)[1] = 2*(Nx_[0]+Nx_[1]) - k;
    }
  else if ( k < (Nx_[0]+1)*(Nx_[1]+1) )
    { // on the interior
      (*i)[0] = (k - numBoundaryPoints_) % (Nx_[0] - 1) + 1;
      (*i)[1] = (k - numBoundaryPoints_) / (Nx_[0] - 1) + 1;
    }
  else
    {
      TEST_FOR_EXCEPTION( true,
	                  std::logic_error,
		          "Illegal running index   k = " << k );
    }

  return i;
}
// =============================================================================
// Defines a series of neighboring boundary nodes.
// This is independent of the actual numbering scheme of the nodes.
Teuchos::RCP<Teuchos::Array<int> >
GridSquare::boundaryPosition ( unsigned int l ) const
{
  int d = 2;
  Teuchos::RCP<Teuchos::Array<int> > i = Teuchos::rcp( new Teuchos::Array<int>(d));

  // start at the bottom left, and go around counter-clockwise
  if (l < Nx_[0])
  { // south
      (*i)[0] = l;
      (*i)[1] = 0;
  }
  else if (l < Nx_[0]+Nx_[1])
  { // east
      (*i)[0] = Nx_[0];
      (*i)[1] = l - Nx_[0];
  }
  else if (l < 2*Nx_[0]+Nx_[1])
  { // north
      (*i)[0] = 2*Nx_[0]+Nx_[1] - l;
      (*i)[1] = Nx_[1];
  }
  else if (l < numBoundaryPoints_)
  { // west
      (*i)[0] = 0;
      (*i)[1] = 2*(Nx_[0]+Nx_[1])  - l;
  }
  else
  {
      TEST_FOR_EXCEPTION( true, std::logic_error,
                          "Given index l=" << l
                          << "larger than the number of boundary nodes n="
                          << numBoundaryPoints_  << ".");
  }
  return i;
}
// =============================================================================
unsigned int
GridSquare::boundaryIndex2globalIndex( unsigned int l ) const
{
	Teuchos::RCP<Teuchos::Array<int> > i = boundaryPosition(l);
	return i2k( i );
}
// =============================================================================
GridSquare::nodeType
GridSquare::getBoundaryNodeType( unsigned int l ) const
{
  Teuchos::RCP<Teuchos::Array<int> > i = boundaryPosition(l);

  if ( (*i)[0]==0 ) {
          if ((*i)[1]==0)
                  return BOTTOMLEFTCONVEX;
          else if ((*i)[1]==(int)Nx_[1])
                  return TOPLEFTCONVEX;
          else
                  return LEFT;
  } else if ( (*i)[0]==(int)Nx_[1] ) {
          if ( (*i)[1]==0 )
                  return BOTTOMRIGHTCONVEX;
          else if ( (*i)[1]==(int)Nx_[1] )
                  return TOPRIGHTCONVEX;
          else
                  return RIGHT;
  } else if ( (*i)[1]==0 )
          return BOTTOM;
  else if ( (*i)[1]==(int)Nx_[1] )
          return TOP;
  else
          return INTERIOR;
}
// =============================================================================
double
GridSquare::cellArea(unsigned int k) const
{
  // TODO Use nodeType association here and possibly move out to Virtual
  if (k == 0 || k == Nx_[0] || k == Nx_[0]+Nx_[1] || k == 2*Nx_[0]+Nx_[1] )
	  // corner
	  return 0.25*h_[0]*h_[1];
  else if (k < numBoundaryPoints_)
	  // edge
      return 0.5*h_[0]*h_[1];
  else
	  // interior
	  return h_[0]*h_[1];
}
// =============================================================================
// maps a 2D index i to a running index k
unsigned int
GridSquare::i2k( Teuchos::RCP<Teuchos::Array<int> > & i) const
{
  int k;

  if ((*i)[1] == 0) // south
      k = (*i)[0];
  else if ((*i)[0] == (int)Nx_[0]) // east
      k = (*i)[1] + Nx_[0];
  else if ((*i)[1] == (int)Nx_[1]) // north
      k = 2*Nx_[0]+Nx_[1] - (*i)[0];
  else if ((*i)[0] == 0) // west
      k = 2*(Nx_[0]+Nx_[1]) - (*i)[1];
  else if ((*i)[0] > 0 && (*i)[0] < (int)Nx_[0] && (*i)[1] > 0 && (*i)[1] < (int)Nx_[1]) // interior
      k = 2*(Nx_[0]+Nx_[1])
        + (Nx_[0] - 1) * ((*i)[1] - 1)
        + (*i)[0] - 1;
  else
      TEST_FOR_EXCEPTION( true, std::logic_error,
	                  "Illegal 2D index   i = " << *i );

  return k;
}
// =============================================================================
// TODO implement this using import/export mechanisms
// TODO templatetize this
Teuchos::RCP<DoubleMultiVector>
GridSquare::permuteGrid2Lexicographic( const DoubleMultiVector & x
		                     ) const
{
  TEST_FOR_EXCEPTION( x.getGlobalLength() != numGridPoints_,
		      std::logic_error,
		      "Global length of the input vector x ("
		      << x.getGlobalLength() << ") does not coincide "
		      << "with with number of unknowns on the grid ("
		     << numGridPoints_ << ")." );

  unsigned int numVectors = x.getNumVectors();
  Teuchos::RCP<DoubleMultiVector> xLexicographic =
                  Teuchos::rcp( new DoubleMultiVector(x.getMap(),numVectors) );

  // Loop through the lexicographic ordering.
  for ( unsigned int l=0; l<numVectors; l++) {
      Teuchos::ArrayRCP<const double> xView = x.getVector(l)->get1dView();
      int k = 0;
      Teuchos::RCP<Teuchos::Array<int> > index = Teuchos::rcp( new Teuchos::Array<int>(2) );
      for (unsigned int j = 0; j < Nx_[1] + 1; j++)
      {
          (*index)[1] = j;
          for (unsigned int i = 0; i < Nx_[0] + 1; i++)
          {
              (*index)[0] = i;
              int kGlobal = i2k(index);
              xLexicographic->replaceGlobalValue( k++, l, xView[kGlobal] );
          }
      }
  }

  return xLexicographic;
}
// =============================================================================
// TODO implement this using import/export mechanisms
// TODO templatetize this
Teuchos::RCP<ComplexMultiVector>
GridSquare::permuteGrid2Lexicographic( const ComplexMultiVector & x
		                             ) const
{
  TEST_FOR_EXCEPTION( x.getGlobalLength() != numGridPoints_,
		      std::logic_error,
		      "Global length of the input vector x ("
		      << x.getGlobalLength() << ") does not coincide "
		      << "with with number of unknowns on the grid ("
		     << numGridPoints_ << ")." );

  unsigned int numVectors = x.getNumVectors();
  Teuchos::RCP<ComplexMultiVector> xLexicographic =
                  Teuchos::rcp( new ComplexMultiVector(x.getMap(),numVectors) );

  // Loop through the lexicographic ordering.
  for ( unsigned int l=0; l<numVectors; l++) {
      Teuchos::ArrayRCP<const std::complex<double> > xView = x.getVector(l)->get1dView();
      int k = 0;
      Teuchos::RCP<Teuchos::Array<int> > index = Teuchos::rcp( new Teuchos::Array<int>(2) );
      for (unsigned int j = 0; j < Nx_[1] + 1; j++)
      {
          (*index)[1] = j;
          for (unsigned int i = 0; i < Nx_[0] + 1; i++)
          {
              (*index)[0] = i;
              int kGlobal = i2k(index);
              xLexicographic->replaceGlobalValue( k++, l, xView[kGlobal] );
          }
      }
  }

  return xLexicographic;
}
// =============================================================================
//// TODO implement this using import/export mechanisms
//// TODO templatetize this
Teuchos::RCP<DoubleMultiVector>
GridSquare::permuteLexicographic2Grid( const DoubleMultiVector & xLexicographic
                                     ) const
{
  TEST_FOR_EXCEPTION( xLexicographic.getGlobalLength() != numGridPoints_,
                      std::logic_error,
                      "Global length of the input vector x ("
                      << xLexicographic.getGlobalLength() << ") does not coincide "
                      << "with with number of unknowns on the grid ("
                      << numGridPoints_ << ")." );

  unsigned int numVectors = xLexicographic.getNumVectors();
  Teuchos::RCP<DoubleMultiVector> x =
                  Teuchos::rcp( new DoubleMultiVector(xLexicographic.getMap(),numVectors) );

  // Loop through the lexicographic ordering.
  for ( unsigned int l=0; l<numVectors; l++) {
      Teuchos::ArrayRCP<const double> xLexicographicView = xLexicographic.getVector(l)->get1dView();
      int k = 0;
      Teuchos::RCP<Teuchos::Array<int> > index = Teuchos::rcp( new Teuchos::Array<int>(2) );
      for (unsigned int j = 0; j < Nx_[1] + 1; j++)
      {
          (*index)[1] = j;
          for (unsigned int i = 0; i < Nx_[0] + 1; i++)
          {
              (*index)[0] = i;
              int kGlobal = i2k(index);
              x->replaceGlobalValue( kGlobal, l, xLexicographicView[k++] );
          }
      }
  }

  return x;
}
// =============================================================================
//// TODO implement this using import/export mechanisms
//// TODO templatetize this
Teuchos::RCP<ComplexMultiVector>
GridSquare::permuteLexicographic2Grid( const ComplexMultiVector & xLexicographic
                                     ) const
{
  TEST_FOR_EXCEPTION( xLexicographic.getGlobalLength() != numGridPoints_,
                      std::logic_error,
                      "Global length of the input vector x ("
                      << xLexicographic.getGlobalLength() << ") does not coincide "
                      << "with with number of unknowns on the grid ("
                      << numGridPoints_ << ")." );

  unsigned int numVectors = xLexicographic.getNumVectors();
  Teuchos::RCP<ComplexMultiVector> x =
                  Teuchos::rcp( new ComplexMultiVector(xLexicographic.getMap(),numVectors) );

  // Loop through the lexicographic ordering.
  for ( unsigned int l=0; l<numVectors; l++) {
      Teuchos::ArrayRCP<const std::complex<double> > xLexicographicView = xLexicographic.getVector(l)->get1dView();
      int k = 0;
      Teuchos::RCP<Teuchos::Array<int> > index = Teuchos::rcp( new Teuchos::Array<int>(2) );
      for (unsigned int j = 0; j < Nx_[1] + 1; j++)
      {
          (*index)[1] = j;
          for (unsigned int i = 0; i < Nx_[0] + 1; i++)
          {
              (*index)[0] = i;
              int kGlobal = i2k(index);
              x->replaceGlobalValue( kGlobal, l, xLexicographicView[k++] );
          }
      }
  }

  return x;
}
// =============================================================================
void
GridSquare::writeWithGrid( const DoubleMultiVector      & x,
                           const Teuchos::ParameterList & params,
                           const std::string            & filePath
                         ) const
{
  Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp(IoFactory::createFileIo(filePath));

  // append grid parameters
  Teuchos::ParameterList extendedParams( params );
  extendedParams.get("scaling", scaling_ );
  extendedParams.get("Nx", Nx_[0] );
  extendedParams.get("Ny", Nx_[1] );

  // reorder the grid to lexicographic ordering
  Teuchos::RCP<DoubleMultiVector> xLexicographic = permuteGrid2Lexicographic(x);

  fileIo->write( *xLexicographic, Nx_, h_, extendedParams);
}
// =============================================================================
void
GridSquare::writeWithGrid( const ComplexMultiVector     & x,
                           const Teuchos::ParameterList & params,
                           const std::string            & filePath
                         ) const
{
  Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp(IoFactory::createFileIo(filePath));

  // append grid parameters
  Teuchos::ParameterList extendedParams( params );
  extendedParams.get("scaling", scaling_ );
  extendedParams.get("Nx", Nx_[0] );
  extendedParams.get("Ny", Nx_[1] );

  // reorder the grid to lexicographic ordering
  Teuchos::RCP<ComplexMultiVector> xLexicographic = permuteGrid2Lexicographic(x);

  fileIo->write( *xLexicographic, Nx_, h_, extendedParams);
}
// =============================================================================
void
GridSquare::read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
                  const std::string                             & filePath,
                  Teuchos::RCP<DoubleMultiVector>               & x,
                  Teuchos::ParameterList                        & params
                )
{
  Teuchos::RCP<IoVirtual> fileIo = Teuchos::RCP<IoVirtual>(
                    IoFactory::createFileIo(filePath));

  Teuchos::RCP<DoubleMultiVector> xLexicographic;
  fileIo->read(Comm, xLexicographic, params);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // create the grid with the just attained information
  TEST_FOR_EXCEPTION( !params.isParameter("Nx"),
                      std::logic_error,
                      "Parameter \"Nx\" not found." );
  TEST_FOR_EXCEPTION( !params.isParameter("Ny"),
                      std::logic_error,
                      "Parameter \"Ny\" not found." );
  Nx_ = Teuchos::tuple<unsigned int>( params.get<unsigned int>("Nx"), params.get<unsigned int>("Ny") );


  TEST_FOR_EXCEPTION( !params.isParameter("scaling"),
                      std::logic_error,
                      "Parameter \"scaling\" not found." );
  scaling_ = params.get<double>("scaling");

  // initialization of the dependent members
  h_                 = Teuchos::tuple<double>( scaling_/Nx_[0], scaling_/Nx_[1] );
  gridDomainArea_    = pow( scaling_, 2 );
  numGridPoints_     = (Nx_[0]+1)*(Nx_[1]+1);
  numBoundaryPoints_ = 2*(Nx_[0]+Nx_[1]);

  // apply the grid ordering
  x = permuteLexicographic2Grid( *xLexicographic );
}
// =============================================================================
