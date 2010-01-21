#include "ioVtk.h"

#include <boost/algorithm/string.hpp>

#include <boost/filesystem/fstream.hpp>


typedef Tpetra::Vector<double,Thyra::Ordinal> DoubleVector;

// =============================================================================
// Constructor
IoVtk::IoVtk(std::string fname) :
  IoVirtual(fname)
{
}
// =============================================================================
// Destructor
IoVtk::~IoVtk()
{
}
// =============================================================================
void
IoVtk::read(const Teuchos::RCP<const Teuchos::Comm<int> > & tComm,
                  Teuchos::RCP<DoubleMultiVector>         & x,
                  Teuchos::ParameterList                  & problemParams) const
{
  // check if file exists
  TEST_FOR_EXCEPTION( !boost::filesystem::exists(fileName_),
	                  std::runtime_error,
			          "File \"" << fileName_ << "\" not found." );

  boost::filesystem::ifstream iFile(fileName_);

  // Don't include ifstream::eofbit and ifstream::failbit as otherwise,
  // getline will throw an exception
  // at the end of the file, while it is actually expected to reach the end of
  // the file.
  iFile.exceptions(std::ifstream::badbit);

  // read the header (in particular the size the of the vectors, point_data)
  int vecSize;
  readVtkHeader(iFile, &vecSize, &problemParams );

  // ---------------------------------------------------------------------------
  // read the vector values

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  TEST_FOR_EXCEPTION( !x.is_null(),
      std::logic_error,
      "Input argument PSI must not point to anything significant.\n"
      << "This is to make sure that the read function can set the map\n"
      << "of the MultiVector.\n"
      << "The error message will be gone as soon as replaceMap is\n"
      << "reimplemented in Trilinos." );

  // create map
  Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > myMap = Teuchos::rcp(new Tpetra::Map<Thyra::Ordinal>(vecSize, 0, tComm));
  int numVectors = 2; // TODO: remove this random value!!!
  x = Teuchos::rcp( new DoubleMultiVector(myMap, numVectors, true));
  ReadScalarsFromVtkFile(iFile, vecSize, x);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // TODO Resurrect the following as soon as replaceMap is back in Trilinos
  // build psi of the entries that we got
  //   if ( psi->getGlobalLength() !=  (unsigned int)NumGlobalElements ) {
  //       // discard all old values, define a new map and plug it in
  //       Teuchos::RCP<Tpetra::Map<int> > newMap
  //           = Teuchos::rcp( new Tpetra::Map<int>( NumGlobalElements, 0, psi->getMap()->getComm() ) );
  //       psi->replaceMap( newMap );
  //   }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  //  // TODO Remove this bit of code and replace by a proper parallel version.
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
IoVtk::read(const Teuchos::RCP<const Teuchos::Comm<int> > & tComm,
                  Teuchos::RCP<ComplexMultiVector>        & x,
                  Teuchos::ParameterList                  & problemParams) const
{
  // check if file exists
  TEST_FOR_EXCEPTION( !boost::filesystem::exists(fileName_),
	                  std::runtime_error,
			          "File \"" << fileName_ << "\" not found." );

  boost::filesystem::ifstream iFile(fileName_);

  // Don't include ifstream::eofbit and ifstream::failbit as otherwise,
  // getline will throw an exception
  // at the end of the file, while it is actually expected to reach the end of
  // the file.
  iFile.exceptions(std::ifstream::badbit);

  // read the parameters
//  ReadParamsFromVtkFile(iFile, problemParams);

  // read the header (in particular the size the of the vectors, point_data)
  int vecSize;
  readVtkHeader( iFile, &vecSize, &problemParams );

  // ---------------------------------------------------------------------------
  // read the vector values
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  TEST_FOR_EXCEPTION( !x.is_null(),
      std::logic_error,
      "Input argument PSI must not point to anything significant.\n"
      << "This is to make sure that the read function can set the map\n"
      << "of the MultiVector.\n"
      << "The error message will be gone as soon as replaceMap is\n"
      << "reimplemented in Trilinos." );

  // create map
  Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > myMap = Teuchos::rcp(new Tpetra::Map<Thyra::Ordinal>(vecSize, 0, tComm));
  int numVectors = 1; // TODO: remove this random value!!!
  x = Teuchos::rcp( new ComplexMultiVector(myMap, numVectors, true));
  ReadScalarsFromVtkFile(iFile, vecSize, x);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Closing
  iFile.close();
}
// =============================================================================
void
IoVtk::write( const DoubleMultiVector              & x,
              const Teuchos::Tuple<unsigned int,2> & Nx,
              const Teuchos::Tuple<double,2>       & h,
              const Teuchos::ParameterList         & problemParams )
{
  boost::filesystem::ofstream vtkfile(fileName_);

  // set the output format
  vtkfile.setf(std::ios::scientific);
  vtkfile.precision(15);

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n";

  // write the parameter list
  writeParameterList(problemParams, vtkfile);

  // write the VTK header
  int numScalars = x.getGlobalLength();
  writeVtkStructuredPointsHeader(vtkfile, Nx, h, numScalars);

  // write the hard data
  writeScalars(x, Nx, vtkfile);

  // close the file
  vtkfile.close();
}
// =============================================================================
void
IoVtk::write( const ComplexMultiVector             & x,
              const Teuchos::Tuple<unsigned int,2> & Nx,
              const Teuchos::Tuple<double,2>       & h,
              const Teuchos::ParameterList         & problemParams )
{
  boost::filesystem::ofstream vtkfile(fileName_);

  // set the output format
  vtkfile.setf(std::ios::scientific);
  vtkfile.precision(15);

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n";

  // write the parameter list
  writeParameterList(problemParams, vtkfile);

  // write the VTK header
  int numScalars = x.getGlobalLength();
  writeVtkStructuredPointsHeader(vtkfile, Nx, h, numScalars);

  // write the hard data
  writeScalars(x, Nx, vtkfile);

  // close the file
  vtkfile.close();
}
// =============================================================================
void
IoVtk::write( const DoubleMultiVector              & x,
              const Teuchos::Tuple<unsigned int,2> & Nx,
              const Teuchos::Tuple<double,2>       & h
            )
{
  // open the file
  boost::filesystem::ofstream vtkfile(fileName_);

  // set the output format
  vtkfile.setf(std::ios::scientific);
  vtkfile.precision(15);

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n";

  // write dummy parameter list
  vtkfile << "# # # # # # # # # # # # # #\n";

  // write the VTK header
  int numScalars = x.getGlobalLength();
  writeVtkStructuredPointsHeader(vtkfile, Nx, h, numScalars);

  // write the hard data
  writeScalars(x, Nx, vtkfile);

  // close the file
  vtkfile.close();
}
// =============================================================================
void
IoVtk::write( const ComplexMultiVector             & x,
              const Teuchos::Tuple<unsigned int,2> & Nx,
              const Teuchos::Tuple<double,2>       & h
            )
{
  // open the file
  boost::filesystem::ofstream vtkfile(fileName_);

  // set the output format
  vtkfile.setf(std::ios::scientific);
  vtkfile.precision(15);

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n";

  // write dummy parameter list
  vtkfile << "# # # # # # # # # # # # # #\n";

  // write the VTK header
  int numScalars = x.getGlobalLength();
  writeVtkStructuredPointsHeader(vtkfile, Nx, h, numScalars);

  // write the hard data
  writeScalars(x, Nx, vtkfile);

  // close the file
  vtkfile.close();
}
// =============================================================================
std::string
IoVtk::strJoin(const std::vector<std::string> & vec, const std::string & sep) const
{
  if (vec.size() == 0)
    return "";

  std::string::size_type size = sep.length() * vec.size();
  for (unsigned int i = 0; i < vec.size(); i++)
    size += vec[i].size();

  std::string tmp;
  tmp.reserve(size);
  tmp = vec[0];
  for (unsigned int i = 1; i < vec.size(); i++)
    tmp = tmp + sep + vec[i];

  return tmp;
}
// =============================================================================
void
IoVtk::writeParameterList(const Teuchos::ParameterList & pList,
                                std::ofstream & ioStream) const
{
  // TODO: Look into dynamically *appending* things to paramStringList to avoid
  // the need for precomputing numEntries.

  // count the number of list entries
  int numEntries = 0;
  Teuchos::map<std::string, Teuchos::ParameterEntry>::const_iterator i;
  for (i = pList.begin(); i != pList.end(); ++i)
    numEntries++;

  ioStream << "PARAMETERS ";

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write the list of parameter values
  std::vector<std::string> paramStringList(numEntries);
  int k = 0;
  for (i = pList.begin(); i != pList.end(); ++i)
    {
	  if (i!=pList.begin() ) {
		  ioStream << ", ";
	  }
      std::string paramName = pList.name(i);
      if (pList.isType<int> (paramName))
    	ioStream <<  "int " << pList.name(i)
    	         << "="
                 << pList.get<int>(paramName);
      else if (pList.isType<unsigned int> (paramName))
      	ioStream <<  "unsigned int " << pList.name(i)
      	         << "="
                 << pList.get<unsigned int>(paramName);
      else if (pList.isType<double> (paramName)){
	char buffer[30];
	sprintf(buffer,"%1.10e",pList.get<double> (paramName));
          ioStream <<  "double " << pList.name(i)
        	       << "="
                   << pList.get<double>(paramName);
      }
      else
        {
          TEST_FOR_EXCEPTION( true,
                              std::logic_error,
                              "Parameter \"" << paramName << "\"is neither of type \"int\", \"unsigned int\", nor \"double\"." );
        }
      k++;
    }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ioStream << " END\n";
}
// =============================================================================
void
IoVtk::writeVtkStructuredPointsHeader( std::ofstream & ioStream,
                                       const Teuchos::Tuple<unsigned int,2>  & Nx,
                                       const Teuchos::Tuple<double,2>        & h,
                                       const int numScalars
                                     ) const
{
  ioStream << "ASCII\n"
		   << "DATASET STRUCTURED_POINTS\n"
		   << "DIMENSIONS " << Nx[0] + 1 << " " << Nx[1] + 1 << " " << 1 << "\n"
		   << "ORIGIN 0 0 0\n"
           << "SPACING " << h[0] << " " << h[1] << " " << 0.0 << "\n"
           << "POINT_DATA " << numScalars << "\n";
}
// =============================================================================
// Note that, when writing the data, the values of psi are assumed to be
// given in lexicographic ordering.
void
IoVtk::writeScalars( const DoubleMultiVector & x,
		             const Teuchos::Tuple<unsigned int,2>  & Nx,
                     std::ofstream & oStream) const
{
  // below this threshold, values are actually printed as 0.0
  double threshold = 1.0e-25;

  int numVectors = x.getNumVectors();
  for (int k = 0; k < numVectors; k++)
    {
      oStream << "SCALARS x" << k << " double\n"
              << "LOOKUP_TABLE default\n";
      Teuchos::ArrayRCP<const double> xKView = x.getVector(k)->get1dView();
      int l = 0;
      for (unsigned int j = 0; j < Nx[1] + 1; j++)
        {
          for (unsigned int i = 0; i < Nx[0] + 1; i++)
            {
              // TODO: Remove this.
              // The following ugly construction makes sure that values as 1.234e-46
              // are actually returned as 0.0. This is necessary as ParaView has
              // issues reading the previous.
              // TODO: Handle this in a more generic fashion.
              double val = xKView[l++];

              double alpha = (fabs(val)<threshold) ? 0.0 : val;
              oStream << alpha << "\n";
            }
        }
    }

}
// =============================================================================
// Note that, when writing the data, the values of psi are assumed to be
// given in lexicographic ordering.
void
IoVtk::writeScalars( const ComplexMultiVector              & z,
		             const Teuchos::Tuple<unsigned int,2>  & Nx,
                           std::ofstream                   & oStream
                   ) const
{
  // below this threshold, values are actually printed as 0.0
  double threshold = 1.0e-25;

  int numVectors = z.getNumVectors();

  for (int k = 0; k < numVectors; k++)
    {
      oStream << "SCALARS psi";
      if (numVectors>1) // add a numbering scheme to psi if necessary
    	  oStream << k;
      oStream << " double 2\n"
              << "LOOKUP_TABLE default\n";
      Teuchos::ArrayRCP<const std::complex<double> > zKView = z.getVector(k)->get1dView();
      int l = 0;
      for (unsigned int j = 0; j < Nx[1] + 1; j++)
        {
          for (unsigned int i = 0; i < Nx[0] + 1; i++)
            {
              // TODO: Remove this.
              // The following ugly construction makes sure that values as 1.234e-46
              // are actually returned as 0.0. This is necessary as ParaView has
              // issues reading the previous.
              // TODO: Handle this in a more generic fashion.
              std::complex<double> val = zKView[l++];

              double re = (fabs(real(val))<threshold) ? 0.0 : real(val);
              oStream << re << " ";

              double im = (fabs(imag(val))<threshold) ? 0.0 : imag(val);
              oStream << im << "\n";
            }
        }
    }

}
// ============================================================================
void
IoVtk::readVtkHeader( std::ifstream          & iFile,
		              int                    * pointData,
		              Teuchos::ParameterList * paramList
		            ) const
{
  std::string buf;

  getline(iFile, buf);
  TEST_FOR_EXCEPTION( buf.compare("# vtk DataFile Version 2.0")!=0,
                      std::logic_error,
                      "Illegal format specification \"" << buf << "\" in file \"" << fileName_ << "\"." );

  // read the parameters
  getline(iFile, buf);
  ReadParamsFromVtkFile(buf, *paramList);


  getline(iFile, buf);
  TEST_FOR_EXCEPTION( buf.compare("ASCII")!=0,
                      std::logic_error,
                      "Illegal format specification \"" << buf << "\" in file \"" << fileName_ << "\"." );

  iFile >> buf;
  TEST_FOR_EXCEPTION( buf.compare("DATASET")!=0,
                      std::logic_error,
                      "Keyword \"DATASET\" not found in file '" << fileName_ << "'."
                      << " Instead found \"" << buf << "\".");

  iFile >> buf;
  TEST_FOR_EXCEPTION( buf.compare("STRUCTURED_POINTS")!=0,
                      std::logic_error,
                      "Keyword \"STRUCTURED_POINTS\" not found in file '" << fileName_ << "'."
                      << " Instead found \"" << buf << "\".");

  iFile >> buf;
  TEST_FOR_EXCEPTION( buf.compare("DIMENSIONS")!=0,
                      std::logic_error,
                      "Keyword \"DIMENSIONS\" not found in file '" << fileName_ << "'."
                      << " Instead found \"" << buf << "\".");
  int n; // read and discard
  iFile >> n;
  iFile >> n;
  iFile >> n;

  iFile >> buf;
  TEST_FOR_EXCEPTION( buf.compare("ORIGIN")!=0,
                      std::logic_error,
                      "Keyword \"ORIGIN\" not found in file '" << fileName_ << "'."
                      << " Instead found \"" << buf << "\".");
  double orig; // read and discard
  iFile >> orig;
  iFile >> orig;
  iFile >> orig;

  iFile >> buf;
  TEST_FOR_EXCEPTION( buf.compare("SPACING")!=0,
                      std::logic_error,
                      "Keyword \"SPACING\" not found in file '" << fileName_ << "'."
                      << " Instead found \"" << buf << "\".");
  double h; // read and discard
  iFile >> h;
  iFile >> h;
  iFile >> h;

  iFile >> buf;
  TEST_FOR_EXCEPTION( buf.compare("POINT_DATA")!=0,
                      std::logic_error,
                      "Keyword \"POINT_DATA\" not found in file '" << fileName_ << "'."
                      << " Instead found \"" << buf << "\".");
  iFile >> *pointData;

  return;
}
// =============================================================================
void
IoVtk::ReadParamsFromVtkFile( const std::string            & paramLine,
                                    Teuchos::ParameterList & fileParams
                            ) const
{

  // buffer the string to be able to *work on it
  std::string aString = paramLine;

  // Assume that the second line in the code has the form
  // "PARAMETERS type1 name1=value1 , type2 name2=value2 , ... END".
  // Iterate over the line and extract what we need:

  // Make sure the line starts with PARAMETERS and ends with END, and erase
  // the two; throw an error if either of the two doesn't exist.
  TEST_FOR_EXCEPTION( !boost::algorithm::starts_with(aString,"PARAMETERS"),
                      std::runtime_error,
                      "Keyword \"PARAMETERS\" missing in file '" << fileName_ << "'" );

  // delete the PARAMETERS string
  aString.erase(0, 10);

  TEST_FOR_EXCEPTION( !boost::algorithm::ends_with(aString,"END"),
                      std::runtime_error,
                      "Keyword \"END\" missing in file '" << fileName_ << "'" );

  // delete END string
  aString.erase(aString.size() - 3, aString.size());

  // trim the string first to avoid dealing with entirely empty parameters lists
  // below
  boost::trim(aString);

  // now, 'explode' the string at the separator sign, ',' (comma)
  std::vector<std::string> strVec;
  boost::algorithm::split(strVec, aString, boost::algorithm::is_any_of(","));

  // Now loop over the separate elements and decode them into the parameter
  // vector.
  Teuchos::Tuple<std::string,3> validDataTypes = Teuchos::tuple<std::string>("int", "unsigned int", "double" );

  for (unsigned int k = 0; k < strVec.size(); k++)
    {

      // trim it!
      boost::trim(strVec[k]);

      // discard empty items
      if (strVec[k].empty())
        break;

      // get the variable type
      int it1 = strVec[k].find(" ");
      std::string type = strVec[k].substr(0, it1);

      bool isValidDataType = false;
      int  typeIndex = -1;
      for ( int l=0; l<validDataTypes.size(); l++ ) {
        it1 = validDataTypes[l].length();
        if ( !strVec[k].substr(0, it1).compare( validDataTypes[l] ) ) {
            typeIndex = l;
            isValidDataType = true;
            break;
        }
      }

      TEST_FOR_EXCEPTION( !isValidDataType,
                          std::runtime_error,
                          "Type \"" << type << "\" is none of the valid data types "
                          << validDataTypes << " in file \"" << fileName_ << "\"." );

      std::string equation = strVec[k].substr(it1);

      // split the rest up by '='
      std::vector<std::string> terms;
      boost::algorithm::split(terms, equation, boost::algorithm::is_any_of("="));

      // trim the variable name
      boost::trim(terms[0]);

      TEST_FOR_EXCEPTION( terms.size()!=2,
          std::logic_error,
          "The expression \"" << equation
          << "\" is not of the form \"psi=5.34\""
          << "in file '" << fileName_ << "'." );

      switch (typeIndex) {
      case(0): // int
      {
          int val = strtol(terms[1].c_str(), NULL, 10);
          fileParams.set(terms[0], val);
          break;
      }
      case(1): // unsigned int
      {
          int val = strtoul(terms[1].c_str(), NULL, 10);
          fileParams.set(terms[0], val);
          break;
      }
      case(2): // double
      {
          double val = strtod(terms[1].c_str(), NULL);
          fileParams.set(terms[0], val);
          break;
      }
      default:
        TEST_FOR_EXCEPTION( true,
            std::logic_error,
            "Type \"" << type << "\" is none of the valid data types "
            << validDataTypes << " in file \"" << fileName_ << "\"." );
      }

    } // loop over parameter list
  return;
}
// ============================================================================
void
IoVtk::ReadScalarsFromVtkFile( std::ifstream                   & iFile,
                               const unsigned int                pointData,
                               Teuchos::RCP<DoubleMultiVector> & scalars
                             ) const
{
  TEUCHOS_ASSERT_EQUALITY( scalars->getGlobalLength(), pointData );
  int numVectors = scalars->getNumVectors();

  // TODO: make the reader work on multiproc environments

  // loop over the SCALARS section and read the thing
  int vecIndex = 0;
  std::string buf;
  iFile >> buf; // read ahead
  while (iFile.good())
    {
	  // read the header
	  int numComponents;
	  readScalarFieldHeader( iFile, buf, numComponents );

      // Now the numbers start rolling.
      // Check that the multivector still has enough room for it.
      TEST_FOR_EXCEPTION( vecIndex+numComponents>numVectors,
                          std::logic_error,
                          "Number of data sets in the VTK larger"
                          << " than what can be stored in the MultiVector "
                          << "(#vectors=" << numVectors << ") in file '" << fileName_ << "'." );

      // Read the elements row by row
      double dummy;
      for (unsigned int k = 0; k < pointData; k++)
        {
          for (int component = 0; component < numComponents; component++)
            {
              iFile >> dummy;
              scalars->replaceLocalValue(k, vecIndex + component, dummy);
            }
        }
      vecIndex += numComponents;

	  iFile >> buf; // read-ahead
	}

  return;
}
// ============================================================================
void
IoVtk::ReadScalarsFromVtkFile( std::ifstream                    & iFile,
 		                       const unsigned int                 pointData,
                               Teuchos::RCP<ComplexMultiVector> & scalars
                             ) const
{
  TEUCHOS_ASSERT_EQUALITY( scalars->getGlobalLength(), pointData );

  int numVectors = scalars->getNumVectors();

  // TODO: make the reader work on multiproc environments

  // Part of the legacy code. See below.
  bool realPart = true;
  Teuchos::RCP<DoubleVector> legacyRealValues = Teuchos::null;

  // loop over the SCALARS section and read the thing
  int vecIndex = 0;
  std::string buf;
  iFile >> buf; // read-ahead

  while (iFile.good())
    {
	  // read the header
	  int numComponents;
	  readScalarFieldHeader( iFile, buf, numComponents );

      // now the numbers start running

  	  if (numComponents==2) {
          // Read the elements row by row
          double re, im;
          for (unsigned int k = 0; k < pointData; k++)
            {
               iFile >> re; iFile >> im;
               scalars->replaceGlobalValue(k, pointData, std::complex<double>(re,im) );
            }
          vecIndex++;
  	  } else if (numComponents==1) { // there's only one component, so assume that the other one is in the next scalar field
  		  // TODO Remove legacy mode.
  		  // LEGACY MODE
  		  // ===========
  		  //
  		  // In former VTK write's the complex data was stores in two SCALARS fields, the first of which one containing
  		  // abs(psi)^2, the second arg(psi). If onle one component is found, assume we have one of those
  		  // legacy files.
  		  //
  	      double dummy;
  	      std::complex<double> alpha;
  	      if (realPart) {
  	  		  std::cout << " ******************************************************************\n"
  					    << " ***\n"
  	  				    << " *** Reading VTK data in legacy mode, suited for GL VTK files\n"
  	  				    << " *** produced before January 2010.\n"
  	  				    << " ***\n"
  	  				    << " ******************************************************************"
  	  				    << std::endl;
  	    	  legacyRealValues = Teuchos::rcp( new DoubleVector(scalars->getMap(),pointData) );
  	  	      for (unsigned int k = 0; k < pointData; k++)
  	  	        {
  	  	           iFile >> dummy;
  	  	           legacyRealValues->replaceGlobalValue(k, dummy );
  	  	        }
  	      } else { // imaginary part
  	    	  // get view of the real part
  	    	  Teuchos::ArrayRCP<const double> legacyRealView = legacyRealValues->get1dView();
	  	      std::complex<double> alpha;
  	  	      for (unsigned int k = 0; k < pointData; k++)
  	  	        {
  	  	           iFile >> dummy;
  	  	           // calculate the complex data from the legacy structure
  	  	           alpha = std::polar( sqrt(legacyRealView[k]), dummy );
  	  	           scalars->sumIntoGlobalValue(k, vecIndex, alpha );
  	  	        }
    	    	vecIndex++; // move on to the next vector -- actually useless as nothing more is expected to come
  	      }
  	      realPart = !realPart; // alternate real and imaginary part
  	  } else {
  		  TEST_FOR_EXCEPTION( true,
  				              std::runtime_error,
  				              "Don't know how to handle three or more components for complex value interpretation." );
  	  }
  	  iFile >> buf; // read-ahead
    }

  return;
}
//==============================================================================
void
IoVtk::readScalarFieldHeader( std::ifstream & iFile,
		                      std::string   & buf,
		                      int           & numComponents
		                    ) const
{
  TEST_FOR_EXCEPTION( buf.compare("SCALARS")!=0,
                      std::logic_error,
                      "Keyword \"SCALARS\" not found in file '" << fileName_ << "'."
                      << " Instead found \"" << buf << "\".");

  iFile >> buf; // name of the variable

  iFile >> buf; // type
  TEST_FOR_EXCEPTION( buf.compare("float")!=0 && buf.compare("double")!=0,
                      std::logic_error,
	                  "None of the keywords \"float\" or \"double\" found in file '" << fileName_ << "'."
	                  << " Instead found \"" << buf << "\".");

  // now follows either nothing, or the number of components
  numComponents = 1; // default if nothing
  getline( iFile, buf ); // get the remainder of the line
  if ( !buf.empty() ) // something here! must be numComponents
      numComponents = atoi( buf.c_str() );

  // make sure the number of components is between 1 and 4
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( numComponents, 1, 5 );

  // LOOKUP_TABLE
  iFile >> buf;
  TEST_FOR_EXCEPTION( buf.compare("LOOKUP_TABLE")!=0,
	                  std::logic_error,
	                  "Keyword \"LOOKUP_TABLE\" not found in file '" << fileName_ << "'."
	                  << " Instead found \"" << buf << "\".");

  // LOOKUP_TABLE name
  iFile >> buf; // discard

}
//==============================================================================
