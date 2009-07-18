class PsiGrid
{
  public:
     PsiGrid( int );
     ~PsiGrid();
     int* k2i( int  );
     int  i2k( int* );

//   protected:
//      // Attributes visible to descendents
  private:
     int Nx;
     int d;
};