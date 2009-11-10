/*
 * Grid.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schlšmer
 */

#include "Grid.h"

#include "glException.h"

#include <EpetraExt_Utils.h> // for toString

// =============================================================================
// Class constructor
Grid::Grid(int nx, double edgeLength) :
  nx_(nx),
  edgeLength_(edgeLength),
  h_(edgeLength / nx)
{
}
// =============================================================================
// Destructor
Grid::~Grid()
{
}
// =============================================================================
int
Grid::getNx() const
{
  return nx_;
}
// =============================================================================
double
Grid::getGridDomainArea() const
{
  return edgeLength_*edgeLength_;
}
// =============================================================================
double
Grid::getEdgeLength() const
{
  return edgeLength_;
}
// =============================================================================
void
Grid::setEdgeLength( const double edgeLength )
{
	edgeLength_ = edgeLength;
	h_          = edgeLength_ / nx_;
}
// =============================================================================
int
Grid::getNumGridPoints() const
{
  // the number of grid points in psi
  return (nx_ + 1) * (nx_ + 1);
}
// =============================================================================
int
Grid::getNumBoundaryPoints() const
{
  return 4*nx_;
}
// =============================================================================
double
Grid::getH() const
{
  return h_;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
Grid::getX(Teuchos::Array<int> i) const
{
  Teuchos::RCP<Teuchos::Array<double> > x = Teuchos::rcp(new Teuchos::Array<
      double>(2));
  (*x)[0] = i[0] * h_;
  (*x)[1] = i[1] * h_;
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
Grid::getXLeft(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[0] -= 0.5 * h_;
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
Grid::getXRight(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[0] += 0.5 * h_;
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
Grid::getXBelow(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[1] -= 0.5 * h_;
  return x;
}
// =============================================================================
Teuchos::RCP<Teuchos::Array<double> >
Grid::getXAbove(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i(k2i(k));
  Teuchos::RCP<Teuchos::Array<double> >x(getX(*i));
  (*x)[1] += 0.5 * h_;
  return x;
}
// =============================================================================
int
Grid::getKLeft( int k ) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[0] -= 1;
  return i2k(i);
}
// =============================================================================
int
Grid::getKRight(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[0] += 1;
  return i2k(i);
}
// =============================================================================
int
Grid::getKBelow(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[1] -= 1;
  return i2k(i);
}
// =============================================================================
int
Grid::getKAbove(int k) const
{
  Teuchos::RCP<Teuchos::Array<int> > i( k2i(k) );
  (*i)[1] += 1;
  return i2k(i);
}
// =============================================================================
// maps a running index k to a 2D index i
Teuchos::RCP<Teuchos::Array<int> >
Grid::k2i(int k) const
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
      std::string message = "Illegal running index   k = " + EpetraExt::toString(k);
      throw glException("Grid::k2i", message);
    }

  return i;
}
// =============================================================================
// Defines a series of neighboring boundary nodes.
// This is independent of the actual numbering scheme of the nodes.
Teuchos::RCP<Teuchos::Array<int> >
Grid::boundaryPosition ( int l ) {
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
			  std::string message = "Given index l=" + EpetraExt::toString(l)
			                      + "larger than the number of border nodes n="
			                      + EpetraExt::toString(4*nx_);
			  throw glException( "Grid::borderNode", message );
		}
		return i;
}
// =============================================================================
int
Grid::boundaryIndex2globalIndex( int l )
{
	Teuchos::RCP<Teuchos::Array<int> > i = boundaryPosition(l);
	return i2k( i );
}
// =============================================================================
Grid::nodeType
Grid::boundaryNodeType( int l )
{
	Teuchos::RCP<Teuchos::Array<int> > i = boundaryPosition(l);

	if ( (*i)[0]==0 ) {
		if ((*i)[1]==0)
			return Grid::BOTTOMLEFTCONVEX;
		else if ((*i)[1]==nx_)
			return Grid::TOPLEFTCONVEX;
		else
			return Grid::LEFT;
	} else if ( (*i)[0]==nx_ ) {
		if ( (*i)[1]==0 )
			return Grid::BOTTOMRIGHTCONVEX;
		else if ( (*i)[1]==nx_ )
			return Grid::TOPRIGHTCONVEX;
		else
			return Grid::RIGHT;
	} else if ( (*i)[1]==0 )
		return Grid::BOTTOM;
	else if ( (*i)[1]==nx_ )
		return Grid::TOP;
	else
		return Grid::INTERIOR;
}
// =============================================================================
double
Grid::cellArea(int k) const
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
Grid::i2k( Teuchos::RCP<Teuchos::Array<int> > & i) const
{
  int k;

  if ((*i)[1] == 0)
    { // south
      k = (*i)[0];
    }
  else if ((*i)[0] == nx_)
    { // east
      k = (*i)[1] + nx_;
    }
  else if ((*i)[1] == nx_)
    { // north
      k = 3 * nx_ - (*i)[0];
    }
  else if ((*i)[0] == 0)
    { // west
      k = 4 * nx_ - (*i)[1];
    }
  else if ((*i)[0] > 0 && (*i)[0] < nx_ && (*i)[1] > 0 && (*i)[1] < nx_)
    { // interior
      k = 4 * nx_ + (nx_ - 1) * ((*i)[1] - 1) + (*i)[0] - 1;
    }
  else
    {
      std::string message = "Illegal 2D index   i = " + Teuchos::toString(*i);
      throw glException("Grid::i2k", message);
    }

  return k;
}
// =============================================================================
// Returns a vector that defines the reordering from a lexicographic grid to
// the ordering present in this grid
void
Grid::lexicographic2grid(std::vector<int> *p) const
{
  // check if for admissible vector size
  unsigned int numUnknowns = (nx_ + 1) * (nx_ + 1);
  if (p->size() != numUnknowns)
    {
      std::string message = "Size of the input vector p ("
          + EpetraExt::toString(int(p->size())) + ") "
          + "does not coincide with with number of unknowns on "
          + " the grid (" + EpetraExt::toString((nx_ + 1) * (nx_ + 1)) + ").";
      throw glException("Grid::lexicographic2grid", message);
    }

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
Grid::reorderToLexicographic( Tpetra::Vector<std::complex<double> > & x
		                    ) const
{
  // check if for admissible vector size
  unsigned int numUnknowns = getNumGridPoints();

  if (x.getGlobalLength() != numUnknowns)
    {
      std::string message = "Global length of the input vector x ("
          + EpetraExt::toString(int(x.getGlobalLength())) + ") "
          + "does not coincide with with number of unknowns on "
          + " the grid (" + EpetraExt::toString(numUnknowns) + ").";
      throw glException("Grid::reorderToLexicographic", message);
    }

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
Grid::reorderFromLexicographic( Tpetra::Vector<std::complex<double> > & x
		                      ) const
{
  // check if for admissible vector size
  unsigned int numUnknowns = getNumGridPoints();

  if (x.getGlobalLength() != numUnknowns)
    {
      std::string message = "Global length of the input vector x ("
          + EpetraExt::toString(int(x.getGlobalLength())) + ") "
          + "does not coincide with with number of unknowns on "
          + " the grid (" + EpetraExt::toString(numUnknowns) + ").";
      throw glException("Grid::reorderToLexicographic", message);
    }

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
