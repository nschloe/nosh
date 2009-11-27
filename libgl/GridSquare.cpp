/*
 * Grid.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schlšmer
 */

#include "GridSquare.h"

#include "ioVirtual.h"
#include "ioFactory.h"

// =============================================================================
// Class constructor
GridSquare::GridSquare(int nx, double scaling) :
  GridVirtual( scaling,
               scaling / nx,
               pow( scaling, 2 ),
               (nx+1)*(nx+1), 4*nx
             ),
  nx_(nx)
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
  (*x)[0] = i[0] * h_;
  (*x)[1] = i[1] * h_;
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getXLeft(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[0] -= 0.5 * h_;
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getXRight(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[0] += 0.5 * h_;
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getXBelow(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[1] -= 0.5 * h_;
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
GridSquare::getXAbove(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[1] += 0.5 * h_;
  return x;
}
// =============================================================================
int
GridSquare::getKLeft( int k ) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[0] -= 1;
  return i2k(i);
}
// =============================================================================
int
GridSquare::getKRight(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[0] += 1;
  return i2k(i);
}
// =============================================================================
int
GridSquare::getKBelow(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[1] -= 1;
  return i2k(i);
}
// =============================================================================
int
GridSquare::getKAbove(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[1] += 1;
  return i2k(i);
}
// =============================================================================
// maps a running index k to a 2D index i
Teuchos::RCP<Teuchos::Array<int> >
GridSquare::k2i(int k) const
{
  int d = 2;
  Teuchos::RCP<Teuchos::Array<int> > i = Teuchos::rcp(
      new Teuchos::Array<int>(d));

  if (k < nx_)
    { // south
      (*i)[0] = k;
      (*i)[1] = 0;
    }
  else if (k < 2 * nx_)
    { // east
      (*i)[0] = nx_;
      (*i)[1] = k - nx_;
    }
  else if (k < 3 * nx_)
    { // north
      (*i)[0] = 3 * nx_ - k;
      (*i)[1] = nx_;
    }
  else if (k < 4 * nx_)
    { // west
      (*i)[0] = 0;
      (*i)[1] = 4 * nx_ - k;
    }
  else if ( k < (nx_+1)*(nx_+1) )
    { // on the interior
      int numBoundaryNodes = 4 * nx_;
      (*i)[0] = (k - numBoundaryNodes) % (nx_ - 1) + 1;
      (*i)[1] = (k - numBoundaryNodes) / (nx_ - 1) + 1;
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
GridSquare::boundaryPosition ( int l ) const
{
  int d = 2;
  Teuchos::RCP<Teuchos::Array<int> > i = Teuchos::rcp( new Teuchos::Array<int>(d));

  // start at the bottom left, and go around counter-clockwise
  if (l < nx_)
  { // south
      (*i)[0] = l;
      (*i)[1] = 0;
  }
  else if (l < 2 * nx_)
  { // east
      (*i)[0] = nx_;
      (*i)[1] = l - nx_;
  }
  else if (l < 3 * nx_)
  { // north
      (*i)[0] = 3 * nx_ - l;
      (*i)[1] = nx_;
  }
  else if (l < 4 * nx_)
  { // west
      (*i)[0] = 0;
      (*i)[1] = 4 * nx_ - l;
  }
  else
  {
      TEST_FOR_EXCEPTION( true, std::logic_error,
                          "Given index l=" << l
                          << "larger than the number of boundary nodes n="
                          << 4*nx_  << ".");
  }
  return i;
}
// =============================================================================
int
GridSquare::boundaryIndex2globalIndex( int l ) const
{
	Teuchos::RCP<Teuchos::Array<int> > i = boundaryPosition(l);
	return i2k( i );
}
// =============================================================================
GridSquare::nodeType
GridSquare::getBoundaryNodeType( int l ) const
{
  Teuchos::RCP<Teuchos::Array<int> > i = boundaryPosition(l);

  if ( (*i)[0]==0 ) {
          if ((*i)[1]==0)
                  return BOTTOMLEFTCONVEX;
          else if ((*i)[1]==nx_)
                  return TOPLEFTCONVEX;
          else
                  return LEFT;
  } else if ( (*i)[0]==nx_ ) {
          if ( (*i)[1]==0 )
                  return BOTTOMRIGHTCONVEX;
          else if ( (*i)[1]==nx_ )
                  return TOPRIGHTCONVEX;
          else
                  return RIGHT;
  } else if ( (*i)[1]==0 )
          return BOTTOM;
  else if ( (*i)[1]==nx_ )
          return TOP;
  else
          return INTERIOR;
}
// =============================================================================
double
GridSquare::cellArea(int k) const
{
  if (k == 0 || k == nx_ || k == 2 * nx_ || k == 3 * nx_)
	  // corner
	  return 0.25*h_*h_;
  else if (k < 4 * nx_)
	  // edge
      return 0.5*h_*h_;
  else
	  // interior
	  return h_*h_;
}
// =============================================================================
// maps a 2D index i to a running index k
int
GridSquare::i2k( Teuchos::RCP<Teuchos::Array<int> > & i) const
{
  int k;

  if ((*i)[1] == 0) // south
      k = (*i)[0];
  else if ((*i)[0] == nx_) // east
      k = (*i)[1] + nx_;
  else if ((*i)[1] == nx_) // north
      k = 3 * nx_ - (*i)[0];
  else if ((*i)[0] == 0) // west
      k = 4 * nx_ - (*i)[1];
  else if ((*i)[0] > 0 && (*i)[0] < nx_ && (*i)[1] > 0 && (*i)[1] < nx_) // interior
      k = 4 * nx_ + (nx_ - 1) * ((*i)[1] - 1) + (*i)[0] - 1;
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
      for (int j = 0; j < nx_ + 1; j++)
      {
          (*index)[1] = j;
          for (int i = 0; i < nx_ + 1; i++)
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
      for (int j = 0; j < nx_ + 1; j++)
      {
          (*index)[1] = j;
          for (int i = 0; i < nx_ + 1; i++)
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
  extendedParams.get("Nx", nx_ );

  // reorder the grid to lexicographic ordering
  Teuchos::RCP<DoubleMultiVector> xLexicographic = permuteGrid2Lexicographic(x);

  fileIo->write( *xLexicographic, nx_, h_, extendedParams);
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
  nx_ = params.get<int>("Nx");

  TEST_FOR_EXCEPTION( !params.isParameter("scaling"),
                      std::logic_error,
                      "Parameter \"scaling\" not found." );
  scaling_ = params.get<double>("scaling");

  // initialization of the dependent members
  h_                 = scaling_ / nx_;
  gridDomainArea_    = pow( scaling_, 2 );
  numGridPoints_     = pow( nx_+1   , 2 );
  numBoundaryPoints_ = 4*nx_;

  // apply the grid ordering
  x = permuteLexicographic2Grid( *xLexicographic );
}
// =============================================================================
