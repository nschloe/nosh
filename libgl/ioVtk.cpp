#include "ioVtk.h"
#include "glException.h"

#include <IOVtkUtils.H> // ParaCont

#include <boost/algorithm/string.hpp>


#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// =============================================================================
// Constructor
IoVtk::IoVtk( std::string fname ):
  IoVirtual(fname)
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
void IoVtk::read( std::vector<double_complex> *psi,
                  Teuchos::ParameterList      *problemParams )
{
  // call ParaCont for parameters
  ReadParamsFromVtkFile( fileName, *problemParams );

  // ---------------------------------------------------------------------------
  // read the vector values
  // TODO: Get rid of this.
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif if

  // TODO: What's *this*?!? Don't make assumptions on the parameter list.
  //       Get rid of it.
  int Nx = problemParams->get<int>("Nx");

  // construct a dummy multivector
  Epetra_Map          StandardMap( (Nx+1)*(Nx+1), 0, Comm );
  Teuchos::RCP<Epetra_MultiVector> tmp
                         = Teuchos::rcp(new Epetra_MultiVector(StandardMap,2) );

  ReadScalarsFromVtkFile( fileName, tmp );

  // build psi of the entries that we got
  psi->resize( (Nx+1)*(Nx+1) );
  for (int k=0; k<(Nx+1)*(Nx+1); k++)
      (*psi)[k] = std::polar( (*tmp)[0][k], (*tmp)[1][k] );

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

}
// =============================================================================



// =============================================================================
void IoVtk::write( const std::vector<double_complex> &psi,
                   const Teuchos::ParameterList      &problemParams,
                   StaggeredGrid                     &sGrid          )
{
  int           Nx = sGrid.getNx(),
                index[2];
  double        h  = sGrid.getH();
  std::ofstream vtkfile;

  // set the output format
  vtkfile.setf( std::ios::scientific );
  vtkfile.precision(15);

  // open the file
  vtkfile.open( fileName.c_str() );

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n";

  // count the number of entries
  int numEntries = 0;
  Teuchos::map<std::string, Teuchos::ParameterEntry>::const_iterator i;
  for (i = problemParams.begin(); i !=problemParams.end(); ++i)
      numEntries++;

  // create the list of parameter values
  std::vector<std::string> paramStringList(numEntries);
  int k=0;
  for (i = problemParams.begin(); i !=problemParams.end(); ++i) {

    std::string type;
    if ( problemParams.isType<int>( problemParams.name(i) ) )
        type = "int";
    else if ( problemParams.isType<double>( problemParams.name(i) ) )
        type = "double";
    else {
        std::string message = "Parameter is neither of type \"int\" not of type \"double\".";
        throw glException( "IoVtk::write",
                           message );
    }

    paramStringList[k] = type + " "
                       + problemParams.name(i)
                       + "="
                       + toString( problemParams.entry(i) );
    k++;
  }

  // print parameter list; join the entries with a comma
  vtkfile << "PARAMETERS ";
  vtkfile << strJoin( paramStringList, " , " );
  vtkfile << " END\n";

  // print the rest of the VTK header
  vtkfile << "ASCII\n"
          << "DATASET STRUCTURED_POINTS\n"
          << "DIMENSIONS " << Nx+1 << " " << Nx+1 << " " << 1 << "\n"
          << "ORIGIN 0 0 0\n"
          << "SPACING " << h << " " << h << " " << 0.0 << "\n"
          << "POINT_DATA " << sGrid.getNumComplexUnknowns() << "\n";

  // Note that, when writing the data, the values of psi are assumed to be
  // given in lexicographic ordering.

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write abs(psi)
  vtkfile << "SCALARS abs(psi) float\n"
          << "LOOKUP_TABLE default\n";
  for (int j=0; j<Nx+1; j++) {
      index[1] = j;
      for (int i=0; i<Nx+1; i++) {
          index[0] = i;
          k = sGrid.i2k( index );
          vtkfile << abs(psi[k]) << "\n";
      }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write arg(psi)
  vtkfile << "SCALARS arg(psi) float\n"
          << "LOOKUP_TABLE default\n";
  for (int j=0; j<Nx+1; j++) {
      index[1] = j;
      for (int i=0; i<Nx+1; i++) {
          index[0] = i;
          k = sGrid.i2k( index );
          vtkfile << arg(psi[k]) << "\n";
      }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // close the file
  vtkfile.close();
}
// =============================================================================


// =============================================================================
std::string IoVtk::strJoin( const std::vector<std::string> & vec,
                            const std::string              & sep  )
{
        if(vec.size()==0)
                return "";
        std::string::size_type size=sep.length()*vec.size();
        for(unsigned int i=0;i<vec.size();i++)
                size+=vec[i].size();

        string tmp;
        tmp.reserve(size);
        tmp=vec[0];
        for(unsigned int i=1;i<vec.size();i++)
                tmp=tmp+sep+vec[i];

        return tmp;
}
// =============================================================================