#include <iostream>
#include "psiGrid.h"

int main(int argc, char *argv[])
{

  // create grid
  int Nx = 10;
  PsiGrid grid(Nx);

  int *i = grid.k2i(55);

  for (int k=0; k<2; k++)
      std::cout << i[k] << std::endl;

}