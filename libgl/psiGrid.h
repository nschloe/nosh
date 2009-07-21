class PsiGrid
{
  public:
     PsiGrid( int nx, double edgelength);

     ~PsiGrid();

     int* k2i( int  );
     int  i2k( int* );

     int getNx();

  private:
     int Nx;
     double Edgelength;
};