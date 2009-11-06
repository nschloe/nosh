#include "ioVtk.h"
#include "ioException.h"

#include <string>

#include <boost/algorithm/string.hpp>

#include <EpetraExt_Utils.h>

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
void IoVtk::read( const Teuchos::RCP<const Teuchos::Comm<int> >        &tComm,
		                Teuchos::RCP<Tpetra::MultiVector<double,int> > &x,
                        Teuchos::ParameterList                         &problemParams
                ) const
{
  std::ifstream iFile;

  // Don't include ifstream::eofbit and ifstream::failbit as otherwise,
  // getline will throw an exception
  // at the end of the file, while it is actually exptected to reach the end of
  // the file.
  iFile.exceptions( std::ifstream::badbit );

  // Opening
  try {
	iFile.open( fileName_.c_str(), std::ios_base::in );
  }
  catch (std::ifstream::failure const &e) {
	std::string message = "Exception opening/reading file '" + fileName_ + "'. "
	+ "Error message '" + e.what() + "'";
	throw ioException( "IoVtk::read", message );
  }

  // read the parameters
  ReadParamsFromVtkFile( iFile, problemParams );

  // read the header (in particular the size the of the vectors, point_data)
  int vecSize = readVtkHeader( iFile );

  // ---------------------------------------------------------------------------
  // read the vector values

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if ( !x.is_null() ) {
      std::string message = std::string("Input argument PSI must not point to anything significant.\n")
                          + std::string("This is to make sure that the read function can set the map\n")
			  + std::string("of the MultiVector.\n")
			  + std::string("The error message will be gone as soon as replaceMap is\n")
			  + std::string("reimplemented in Trilinos.");
      throw ioException( "IoVtk::read", message );
  }

  // create map
  Teuchos::RCP<const Tpetra::Map<int> > myMap =
		              Teuchos::rcp(new Tpetra::Map<int>( vecSize, 0, tComm ) );
  int numVectors = 2; // TODO: remove this random value!!!
  x = Teuchos::rcp( new Tpetra::MultiVector<double,int>(myMap,numVectors,true) );
  ReadScalarsFromVtkFile( iFile, x );
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


//  // TODO:
//  // Remove this bit of code and replace by a proper parallel version.
//  for (unsigned int k=0; k<psi->getLocalLength(); k++) {
//      int kGlobal = psi->getMap()->getGlobalElement(k);
//      double_complex z = std::polar( (*tmp)[0][kGlobal], (*tmp)[1][kGlobal] );
//      psi->replaceLocalValue( k, z );
//  }


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
//   ReadScalarsFromVtkFile( fileName_, scalars );


  // Closing
  iFile.close();
}
// =============================================================================
void
IoVtk::write( const Tpetra::MultiVector<double,int> & x,
              const int                               Nx,
              const double                            h,
              const Teuchos::ParameterList          & problemParams
            ) const
{
  std::ofstream vtkfile;
  
  // set the output format
  vtkfile.setf( std::ios::scientific );
  vtkfile.precision(15);

  // open the file
  vtkfile.open( fileName_.c_str() );

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n";

  // write the parameter list
  writeParameterList( problemParams, vtkfile );

  // write the VTK header
  int numScalars = x.getGlobalLength();
  writeVtkStructuredPointsHeader( vtkfile, Nx, h, numScalars );

  // write the hard data
  writeScalars( x, Nx, vtkfile );

  // close the file
  vtkfile.close();
}
// =============================================================================
void
IoVtk::write( const Tpetra::MultiVector<double,int> & x,
              const int                               Nx,
              const double                            h
            ) const
{
  std::ofstream vtkfile;
  
  // set the output format
  vtkfile.setf( std::ios::scientific );
  vtkfile.precision(15);

  // open the file
  vtkfile.open( fileName_.c_str() );

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n";

  // write dummy parameter list
  vtkfile << "# # # # # # # # # # # # # #\n";

  // write the VTK header
  int    numScalars = x.getGlobalLength();
  writeVtkStructuredPointsHeader( vtkfile, Nx, h, numScalars );

  // write the hard data
  writeScalars( x, Nx, vtkfile );

  // close the file
  vtkfile.close();
}
// =============================================================================
std::string
IoVtk::strJoin( const std::vector<std::string> & vec,
                const std::string              & sep
              ) const
{
  if(vec.size()==0)
      return "";

  std::string::size_type size=sep.length()*vec.size();
  for(unsigned int i=0;i<vec.size();i++)
      size+=vec[i].size();

  std::string tmp;
  tmp.reserve(size);
  tmp=vec[0];
  for(unsigned int i=1;i<vec.size();i++)
      tmp=tmp+sep+vec[i];

  return tmp;
}
// =============================================================================
void
IoVtk::writeParameterList( const Teuchos::ParameterList & pList,
                           std::ofstream                & ioStream
                         ) const
{
  // TODO: Look into dynamically *appending*things to paramStringList to avoid
  // the need for precomputing numEntries.

  // count the number of list entries
  int numEntries = 0;
  Teuchos::map<std::string, Teuchos::ParameterEntry>::const_iterator i;
  for (i = pList.begin(); i!=pList.end(); ++i)
      numEntries++;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // create the list of parameter values
  std::vector<std::string> paramStringList(numEntries);
  int k=0;
  for (i = pList.begin(); i!=pList.end(); ++i) {
    std::string paramName = pList.name(i);
    if ( pList.isType<int>( paramName ) )
        paramStringList[k] = "int "
                           + pList.name(i)
                           + "="
                           + EpetraExt::toString( pList.get<int>(paramName) );
    else if ( pList.isType<double>( paramName ) )
        paramStringList[k] = "double "
                           + pList.name(i)
                           + "="
                           + EpetraExt::toString( pList.get<double>(paramName) );
    else {
        std::string message =
                 "Parameter is neither of type \"int\" not of type \"double\".";
        throw ioException( "IoVtk::writeParameterList",
                           message );
    }
    k++;
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ioStream << "PARAMETERS ";
  ioStream << strJoin( paramStringList, " , " );
  ioStream << " END\n";
}
// =============================================================================
void
IoVtk::writeVtkStructuredPointsHeader( std::ofstream & ioStream,
                                       const int     Nx,
                                       const double  h,
                                       const int     numScalars
                                     ) const
{
  ioStream << "ASCII\n"
          << "DATASET STRUCTURED_POINTS\n"
          << "DIMENSIONS " << Nx+1 << " " << Nx+1 << " " << 1 << "\n"
          << "ORIGIN 0 0 0\n"
          << "SPACING " << h << " " << h << " " << 0.0 << "\n"
          << "POINT_DATA " << numScalars << "\n";
}
// =============================================================================
void
IoVtk::writeScalars( const Tpetra::MultiVector<double,int> & x,
                     const int                               Nx,
                           std::ofstream                   & oStream
                   ) const
{
  // Note that, when writing the data, the values of psi are assumed to be
  // given in lexicographic ordering.

  int numVectors = x.getNumVectors();
  for (int k=0; k<numVectors; k++ ) {
      oStream << "SCALARS x" << k << " float\n"
              << "LOOKUP_TABLE default\n";
      Teuchos::ArrayRCP<const double> xKView = x.getVector(k)->get1dView();
      int l = 0;
      for (int j=0; j<Nx+1; j++) {
          for (int i=0; i<Nx+1; i++) {
        	  // TODO: Remove this.
              // The following ugly construction makes sure that values as 1.234e-46
              // are actually returned as 0.0. This is necessary as ParaView has
              // issues reading the previous.
              // TODO: Handle this in a more generic fashion.
              double val = xKView[l++];
              if (fabs(val)<1.0e-25)
                  oStream << 0.0 << "\n";
              else
                  oStream << val << "\n";
          }
      }
  }

}
// =============================================================================
bool
IoVtk::ReadParamsFromVtkFile( std::ifstream    &iFile,
                              Teuchos::ParameterList &fileParams
                            ) const
{
   std::string fname   = "ReadParamsFromVtkFile";

   // Dummy string
   std::string aString;

   getline(iFile,aString);
   getline(iFile,aString); // second line

   // Assume that the second line in the code has the form
   // "PARAMETERS type1 name1=value1 , type2 name2=value2 , ... END".
   // Iterate over the line and extract what we need:

   // Make sure the line starts with PARAMETERS and ends with END, and erase
   // the two; throw an error if either of the two doesn't exist.
   if ( boost::algorithm::starts_with(aString,"PARAMETERS") ) {
     aString.erase(0,10);
   } else {
     std::string message = "Keyword \"PARAMETERS\" missing.";
     throw ioException( fname, message );
   }
   if ( boost::algorithm::ends_with(aString,"END") ) {
     aString.erase(aString.size()-3,aString.size());
   } else {
     std::string message = "Keyword \"END\" missing.";
     throw ioException( fname, message );
   }

   // now, 'explode' the string at the separator sign, ',' (comma)
   std::vector<std::string> strVec;
   boost::algorithm::split( strVec,
                            aString,
                            boost::algorithm::is_any_of(",") );

  // Now loop over the separate elements and decode them into the parameter
  // vector.
  for (unsigned int k=0; k<strVec.size(); k++) {

    // trim it!
    boost::trim(strVec[k]);

    // get the variable type
    int it1 = strVec[k].find(" ");
    std::string type = strVec[k].substr(0,it1);

    if ( type.compare("int") &&
      type.compare("double") ) {
      std::string message = "Type \"" + type + "\" is neither \"int\""
      + " nor \"double\".";
    throw ioException( fname, message );
    }

    std::string equation = strVec[k].substr(it1);

    // split the rest up by '='
    std::vector<std::string> terms;
    boost::algorithm::split( terms,
                             equation,
                             boost::algorithm::is_any_of("=") );

    // trim the variable name
    boost::trim( terms[0] );


    if (terms.size()!=2) {
      std::string message = std::string("The expression \"") + equation
                          + std::string("\" is not of the type \"psi=5.34\".");
      throw ioException( fname, message );
    }

    if ( type.compare("int")==0 ) {
      int val = strtol( terms[1].c_str(), NULL, 10 );
      fileParams.set( terms[0], val );
    } else if ( type.compare("double")==0 ) {
      double val = strtod( terms[1].c_str(), NULL );
      fileParams.set( terms[0], val );
    } else {
      std::string message = "Type \"" + type + "\" is neither \"int\""
      + " nor \"double\".";
      throw ioException( fname, message );
    }
  }

  return true;
}
// ============================================================================
int
IoVtk::readVtkHeader( std::ifstream    &iFile
		            ) const
{
  // check format argument
  std::string buf;
  getline(iFile,buf);
  if ( buf.compare("ASCII") ) {
	std::string message = "VTK file not encoded in \"ASCII\".";
	throw ioException( "IoVtk::readVtkHeader", message );
  }

  // get POINT_DATA
  while (buf.find("POINT_DATA")==std::string::npos)
	getline(iFile,buf);
  int startNum = buf.find_first_of("0123456789");
  int pointData = strtol( buf.substr(startNum).c_str(),
						  NULL, 10 );

  return pointData;
}
// ============================================================================
bool
IoVtk::ReadScalarsFromVtkFile( std::ifstream                                  & iFile,
                               Teuchos::RCP<Tpetra::MultiVector<double,int> > & scalars
                             ) const
{
  std::string fname = "IoVtk::ReadScalarsFromVtkFile";

  int numVectors = scalars->getNumVectors();

  // TODO: make the reader work on multiproc environments

  // loop over the SCALARS section and read the thing
  int vecIndex = 0;
  std::string buf;
  while ( !iFile.eof() ) {
    // move forward to the next SCALARS line
    // TODO: come up with something more elegant here
    bool scalarsFound = false;
    // TODO: With ifstream::eofbit or ifstream::failbit, getline will throw
    //       an exception here. Find out why, and fix.
    while (!iFile.eof()) {
      getline(iFile,buf);
      if (buf.find("SCALARS")!=std::string::npos) {
        scalarsFound = true;
        break;
      }
    }
    if (!scalarsFound) break;

    // skip LOOKUP_TABLE line
    getline(iFile,buf);

    // TODO: make sure we're dealing with FLOATS here

    // Read the number of components; if none is given, take the default (1).
    int numComponents = 1;
    unsigned int startNum = buf.find_first_of("012456789");
    if ( startNum!=std::string::npos )
      numComponents = strtol( buf.substr(startNum,buf.size()-startNum).c_str(),
                              NULL, 10 );

    // check that the multivector still has enough room for it
    if ( vecIndex+numComponents>numVectors ) {
      std::string message = std::string("Number of data sets in the VTK larger")
      + " than what can be stored in the MultiVector "
      + "(#vectors=" + EpetraExt::toString(numVectors)
      + ").";
      throw ioException( fname, message );
    }

    // Read the elements row by row
    double dummy;
    int vecSize = scalars->getGlobalLength();
    for ( int k=0; k<vecSize; k++ ){
      for ( int component=0; component<numComponents; component++ ) {
        iFile >> dummy;
        scalars->replaceLocalValue( k, vecIndex+component, dummy );
      }
    }
    vecIndex += numComponents;
  }

  return true;
}
//==============================================================================
