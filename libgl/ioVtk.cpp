#include "ioVtk.h"
#include "glException.h"

#include <IOVtkUtils.H> // ParaCont

#include <boost/algorithm/string.hpp>

#include <Epetra_Map.h>

#include <EpetraExt_Utils.h>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// TODO:
// Remove this:
#include <Teuchos_Comm.hpp>

// =============================================================================
// Constructor
IoVtk::IoVtk( std::string fname ):
  IoVirtual(fname)
{
}
// =============================================================================
// Destructor
IoVtk::~IoVtk()
{
}
// =============================================================================
void IoVtk::read( Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi,
		  Teuchos::RCP<Teuchos::Comm<int> >                      comm,
                  Teuchos::ParameterList                                 *problemParams )
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
#endif

  // TODO: What's *this*?!? Don't make assumptions on the parameter list.
  //       Get rid of it.
  int Nx = problemParams->get<int>("Nx");
  int NumGlobalElements = (Nx+1)*(Nx+1);
  
  // construct a dummy multivector
  Epetra_Map          StandardMap( NumGlobalElements, 0, Comm );
  Teuchos::RCP<Epetra_MultiVector> tmp
                         = Teuchos::rcp(new Epetra_MultiVector(StandardMap,2) );

  ReadScalarsFromVtkFile( fileName, tmp );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if ( !psi.is_null() ) {
      std::string message = std::string("Input argument PSI must not point to anything significant.\n")
                          + std::string("This is to make sure that the read function can set the map\n")
			  + std::string("of the MultiVector.\n")
			  + std::string("The error message will be gone as soon as replaceMap is\n")
			  + std::string("reimplemented in Trilinos.");
      throw glException( "IoVtk::read", message );
  }
  
  // define map
  Teuchos::RCP<Tpetra::Map<int> > newMap
          = Teuchos::rcp( new Tpetra::Map<int>( NumGlobalElements, 0, comm ) );
  psi = Teuchos::rcp( new Tpetra::MultiVector<double_complex,int>(newMap,1) );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // TODO:
  // Resurrect the following as soon as replaceMap is back in Trilinos
  // build psi of the entries that we got
//   if ( psi->getGlobalLength() !=  (unsigned int)NumGlobalElements ) {
//       // discard all old values, define a new map and plug it in
//       Teuchos::RCP<Tpetra::Map<int> > newMap
//           = Teuchos::rcp( new Tpetra::Map<int>( NumGlobalElements, 0, psi->getMap()->getComm() ) );
//       psi->replaceMap( newMap );
//   }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  // TODO:
  // Remove this bit of code and replace by a proper parallel version.
  for (unsigned int k=0; k<psi->getLocalLength(); k++) {
      int kGlobal = psi->getMap()->getGlobalElement(k);
      double_complex z = std::polar( (*tmp)[0][kGlobal], (*tmp)[1][kGlobal] );
      psi->replaceLocalValue( k, 0, z );
  }

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
void IoVtk::write( const Tpetra::MultiVector<double_complex,int> &psi,
                   const Teuchos::ParameterList                  &problemParams,
                   StaggeredGrid                                 &sGrid          )
{
  int           Nx = sGrid.getNx();
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

    std::string paramName = problemParams.name(i);
    if ( problemParams.isType<int>( paramName ) )
        paramStringList[k] = "int "
                           + problemParams.name(i)
                           + "="
                           + EpetraExt::toString( problemParams.get<int>(paramName) );
    else if ( problemParams.isType<double>( paramName ) )
        paramStringList[k] = "double "
                           + problemParams.name(i)
                           + "="
                           + EpetraExt::toString( problemParams.get<double>(paramName) );
    else {
        std::string message = "Parameter is neither of type \"int\" not of type \"double\".";
        throw glException( "IoVtk::write",
                           message );
    }
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

  Teuchos::ArrayRCP<const double_complex> psiView = psi.getVector(0)->get1dView();
  Teuchos::Array<int> index(2);
  for (int j=0; j<Nx+1; j++) {
      index[1] = j;
      for (int i=0; i<Nx+1; i++) {
          index[0] = i;
          k = sGrid.i2k( index );
          // The following ugly construction makes sure that values as 1.234e-46
          // are actually returned as 0.0. This is necessary as ParaView has
          // issues reading the previous.
          // TODO: Handle this in a more generic fashion.
          double val = abs(psiView[k]);
          if (val<1.0e-25) {
              vtkfile << 0.0 << "\n";
          } else {
              vtkfile << abs(psiView[k]) << "\n";
          }
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
          vtkfile << arg(psiView[k]) << "\n";
      }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // close the file
  vtkfile.close();
}
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