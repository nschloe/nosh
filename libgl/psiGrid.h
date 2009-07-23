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

  private:
     int Nx; //!< Number of grid pieces in both x- and y-direction
};