#include "ioVtk.h"

#include <IOVtkUtils.H> // ParaCont

// =============================================================================
// Constructor
IoVtk::IoVtk( StaggeredGrid &sGrid ):
  SGrid(sGrid) // make this &sGrid
{
}
// =============================================================================



// =============================================================================
// Destructor
IoVtk::~IoVtk()
{
}
// =============================================================================


// =============================================================================
void IoVtk::read( const std::string           &fileName,
                  std::vector<double_complex> *psi,
                  Teuchos::ParameterList      *problemParams )
{

  // call ParaCont for parameters
  ReadParamsFromVtkFile( fileName, *problemParams );

  // call ParaCont for scalars

//   // ALL THIS STUFF HAS TO GO ONE psi IS AN EPETRA_VECTOR!!!
//   // create a standard map for one core
//   Epetra_SerialComm Comm;
//   int MyPID = Comm.MyPID();
//   int numMyUnknowns;
//   int Nx = problemParams->get<int>("Nx");
//   if (MyPID==0)
//       numMyUnknowns = (Nx+1)*(Nx+1);
//   else
//       numMyUnknowns = 0;
// 
//   int NumGlobalElements = (Nx+1)*(Nx+1);
//   Epetra_Map* StandardMap = new Epetra_Map( NumGlobalElements, 0, Comm );
// 
//   Teuchos::RCP<Epetra_MultiVector> scalars = Teuchos::RCP(Epetra_MultiVector(StandardMap,1));
//   ReadScalarsFromVtkFile( fileName, scalars );


  // -------------------------------------------------------------------------
  // The following is a more complete approach to reading the parameters than
  // what is implemented in ParaCont, but unfinished. Suspend for now.
  //
  // One major difference with ParaCont's implementation is that this routine
  // does not assume anything about problemParams at input, and specifically
  // doesn't check against non-existing list entries.
  // -------------------------------------------------------------------------
//   std::ifstream iFile;
//   iFile.exceptions (   std::ifstream::eofbit
//                      | std::ifstream::failbit
//                      | std::ifstream::badbit  );
// 
//   // open the file
//   try {
//       iFile.open( fileName.c_str(), std::ios_base::in );
//   }
//   catch (std::ifstream::failure e) {
//     std::cout << "Exception opening file " << fileName << std::endl;
//   }
// 
//   // Dummy string
//   std::string aString;
// 
//   // Opening a file to check whether a valid list of parameters
//   // is given: the file will be opened again subsequently
//   getline(iFile,aString);
//   getline(iFile,aString); // second line
// 
//   // Assume that the second line in the code has *exactly* the form
//   // "PARAMETERS name1=value1 name2=value2 ... END".
//   // Iterate over the line and extract what we need:
//   int it      = 0,
//       innerIt = aString.find(" ");
// 
//   // find the first string and make sure it's "PARAMETERS"
//   if ( aString.substr(it,innerIt).compare("PARAMETERS") ) {
//       std::cerr << "Parameter list does not start with \"PARAMETERS\". Abort."
//                 << std::endl;
//       exit(EXIT_FAILURE);
//   }
//   innerIt++;
// 
//   // scan through the rest of the line
//   while ( innerIt<aString.size() ) {
//       it = innerIt;
//       innerIt = aString.find(" ",it);
//       if ( innerIt==string::npos ) break;// nothing was found
// 
//       // itemString must be of the form "xxxx=341.21" or similar -- parse it
//       std::string itemString = aString.substr(it,innerIt-it);
//       int it1 = itemString.find("=");
//       if ( it1==string::npos ) {
//           std::cerr << "Parameter item not of the form \"psi=5.34\". Abort."
//                     << std::endl;
//           exit(EXIT_FAILURE);
//       }
// std::cout << ">>" << itemString.substr(0,it1) << "<<" << std::endl;
// std::cout << ">>" << itemString.substr(it1+1,itemString.size()-it1) << "<<" << std::endl;
// 
//       innerIt++;
//   }
// 
//   // check if the last word is indeed "END"
//   if ( aString.substr(it,aString.size()).compare("END") ) {
//           std::cerr << "Parameter list does not end with \"END\". Abort."
//                     << std::endl;
//           exit(EXIT_FAILURE);
//   }
// exit(0);
  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

}
// =============================================================================




// =============================================================================
void IoVtk::write( const std::string                 &filename,
                   const std::vector<double_complex> &psi,
                   const Teuchos::ParameterList      &problemParams )
{
  int           Nx = SGrid.getNx(),
                k,
                index[2];
  double        h  = SGrid.getH();
  std::ofstream vtkfile;

  // open the file
  vtkfile.open( filename.c_str() );

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n";

  // print parameter list
  vtkfile << "PARAMETERS ";
  Teuchos::map<std::string, Teuchos::ParameterEntry>::const_iterator i;
  for (i = problemParams.begin(); i !=problemParams.end(); ++i)
    vtkfile << problemParams.name(i) << "=" << problemParams.entry(i) << " ";
  vtkfile << "END\n";

  // print the rest of the VTK header
  vtkfile << "ASCII\n"
          << "DATASET STRUCTURED_POINTS\n"
          << "DIMENSIONS " << Nx+1 << " " << Nx+1 << " " << 1 << "\n"
          << "ORIGIN 0 0 0\n"
          << "SPACING " << h << " " << h << " " << 0.0 << "\n"
          << "POINT_DATA " << SGrid.getNumComplexUnknowns() << "\n";

  // Note that, when writing the data, the values of psi are assumed to be
  // given in lexicographic ordering.

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write abs(psi)
  vtkfile << "SCALARS abs(psi) float\n"
          << "LOOKUP_TABLE default\n";
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          vtkfile << abs(psi[k]) << "\n";
      }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write arg(psi)
  vtkfile << "SCALARS arg(psi) float\n"
          << "LOOKUP_TABLE default\n";
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          vtkfile << arg(psi[k]) << "\n";
      }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // close the file
  vtkfile.close();
}
// =============================================================================