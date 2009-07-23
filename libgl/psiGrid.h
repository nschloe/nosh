/********************************************//**
 * Grid for \f$\psi\f$.
 ***********************************************/
class PsiGrid
{
  public:
     PsiGrid( int nx );

     ~PsiGrid();

     int getNx(); //!< returns Nx

//     int* k2i( int  ); //!< Converts a running index k to a grid index i
     int  i2k( int* ); //!< Converts a grid index i to a running index k

     /*! Indicates whether a node sits in a corner of the domain, on an edge,
         or strictly inside it. */
     enum nodeType { CORNER, EDGE, INTERIOR };

     /*! For a given node number k, indicates whether the node sits in a corner
         of the domain, on an edge, or strictly inside it. */
     nodeType k2nodeType( int k );

  private:
     int Nx; //!< Number of grid pieces in both x- and y-direction
};