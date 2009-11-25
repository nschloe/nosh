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
GridSquare::GridSquare(int nx, double edgeLength) :
  nx_(nx),
  edgeLength_(edgeLength),
  h_(edgeLength / nx)
{
}
// =============================================================================
// Destructor
GridSquare::~GridSquare()
{
}
// =============================================================================
int
GridSquare::getNx() const
{
  return nx_;
}
// =============================================================================
double
GridSquare::getGridDomainArea() const
{
  return edgeLength_*edgeLength_;
}
// =============================================================================
double
GridSquare::getEdgeLength() const
{
  return edgeLength_;
}
// =============================================================================
void
GridSquare::setEdgeLength( const double edgeLength )
{
	edgeLength_ = edgeLength;
	h_          = edgeLength / nx_;
}
// =============================================================================
int
GridSquare::getNumGridPoints() const
{
  // the number of grid points in psi
  return (nx_ + 1) * (nx_ + 1);
}
// =============================================================================
int
GridSquare::getNumBoundaryPoints() const
{
  return 4*nx_;
}
// =============================================================================
double
GridSquare::getH() const
{
  return h_;
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
GridSquare::boundaryPosition ( int l ) {
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
					            << "larger than the number of border nodes n="
					            << 4*nx_ );
		}
		return i;
}
// =============================================================================
int
GridSquare::boundaryIndex2globalIndex( int l )
{
	Teuchos::RCP<Teuchos::Array<int> > i = boundaryPosition(l);
	return i2k( i );
}
// =============================================================================
GridSquare::nodeType
GridSquare::boundaryNodeType( int l )
{
	Teuchos::RCP<Teuchos::Array<int> > i = boundaryPosition(l);

	if ( (*i)[0]==0 ) {
		if ((*i)[1]==0)
			return GridSquare::BOTTOMLEFTCONVEX;
		else if ((*i)[1]==nx_)
			return GridSquare::TOPLEFTCONVEX;
		else
			return GridSquare::LEFT;
	} else if ( (*i)[0]==nx_ ) {
		if ( (*i)[1]==0 )
			return GridSquare::BOTTOMRIGHTCONVEX;
		else if ( (*i)[1]==nx_ )
			return GridSquare::TOPRIGHTCONVEX;
		else
			return GridSquare::RIGHT;
	} else if ( (*i)[1]==0 )
		return GridSquare::BOTTOM;
	else if ( (*i)[1]==nx_ )
		return GridSquare::TOP;
	else
		return GridSquare::INTERIOR;
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
// Returns a vector that defines the reordering from a lexicographic grid to
// the ordering present in this grid
void
GridSquare::lexicographic2grid(std::vector<int> *p) const
{
  // check if for admissible vector size
  unsigned int numUnknowns = (nx_ + 1) * (nx_ + 1);

  TEST_FOR_EXCEPTION( p->size() != numUnknowns,
		              std::logic_error,
		              "Size of the input vector p (" << p->size() <<") "
		              << "does not coincide with with number of unknowns on "
		              << " the grid (" << (nx_+1)*(nx_+1)  << ")." );

  int k = 0;
  Teuchos::RCP<Teuchos::Array<int> > index = Teuchos::rcp( new Teuchos::Array<int>(2) );
  for (int j = 0; j < nx_ + 1; j++)
    {
      (*index)[1] = j;
      for (int i = 0; i < nx_ + 1; i++)
        {
          (*index)[0] = i;
          (*p)[k++] = i2k(index);
        }
    }

}
// =============================================================================
// TODO implement this using import/export mechanisms
// TODO templatetize this
void
GridSquare::reorderToLexicographic( Tpetra::Vector<std::complex<double> > & x
		                    ) const
{
  // check if for admissible vector size
  unsigned int numUnknowns = getNumGridPoints();

  TEST_FOR_EXCEPTION( x.getGlobalLength() != numUnknowns,
		              std::logic_error,
		              "Global length of the input vector x ("
		              << x.getGlobalLength() << ") does not coincide "
		              << "with with number of unknowns on the grid ("
		              << numUnknowns << ")." );

  // Make a temporary copy of the full vector.
  Tpetra::Vector<std::complex<double> > tmp( x );

  // Loop through the lexicographic ordering.
  // This really depends on the grid to be rectangular.
  Teuchos::ArrayRCP<const std::complex<double> > tmpView = tmp.get1dView();
  int k = 0;
  Teuchos::RCP<Teuchos::Array<int> > index = Teuchos::rcp( new Teuchos::Array<int>(2) );
  for (int j = 0; j < nx_ + 1; j++)
    {
      (*index)[1] = j;
      for (int i = 0; i < nx_ + 1; i++)
        {
          (*index)[0] = i;
          int kGlobal = i2k(index);
          x.replaceGlobalValue( kGlobal, tmpView[k++] );
        }
    }

}
// =============================================================================
// TODO implement this using import/export mechanisms
// TODO templatetize this
void
GridSquare::reorderFromLexicographic( Tpetra::Vector<std::complex<double> > & x
		                      ) const
{
  // check if for admissible vector size
  unsigned int numUnknowns = getNumGridPoints();

  TEST_FOR_EXCEPTION( x.getGlobalLength() != numUnknowns,
		              std::logic_error,
		              "Global length of the input vector x ("
		              << x.getGlobalLength() << ") does not coincide "
		              << "with with number of unknowns on the grid ("
		              << numUnknowns << ")." );

  // Make a temporary copy of the full vector.
  Tpetra::Vector<std::complex<double> > tmp( x );

  // Loop through the lexicographic ordering.
  // This really depends on the grid to be rectangular.
  Teuchos::ArrayRCP<const std::complex<double> > tmpView = tmp.get1dView();
  int k = 0;
  Teuchos::RCP<Teuchos::Array<int> > index = Teuchos::rcp( new Teuchos::Array<int>(2) );
  for (int j = 0; j < nx_ + 1; j++)
    {
      (*index)[1] = j;
      for (int i = 0; i < nx_ + 1; i++)
        {
          (*index)[0] = i;
          int kGlobal = i2k(index);
          x.replaceGlobalValue( k++, tmpView[kGlobal] );
        }
    }

}
// =============================================================================
void
GridSquare::writeWithGrid( const Tpetra::MultiVector<double,int> & x,
                     const Teuchos::ParameterList &params,
                     const std::string &filePath) const
{
  Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp(IoFactory::createFileIo(filePath));
  fileIo->write( x, nx_, h_, params);
}
// =============================================================================
// ATTENTION: Not a member of GridSquare!
void
readWithGrid( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
              const std::string                             & filePath,
              Teuchos::RCP<Tpetra::MultiVector<double> >    & x,
              Teuchos::RCP<GridSquare>                      & grid,
              Teuchos::ParameterList                        & params )
{
  Teuchos::RCP<IoVirtual> fileIo = Teuchos::RCP<IoVirtual>(
                    IoFactory::createFileIo(filePath));
  fileIo->read(Comm, x, params);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // create the grid with the just attained information
  int    Nx         = params.get<int>("Nx");
  double edgeLength = params.get<double>("edge length");
  grid = Teuchos::rcp( new GridSquare( Nx, edgeLength) );
}
// =============================================================================
