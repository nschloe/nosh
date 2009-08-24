#include "ioVirtual.h"

class IoVtk: public IoVirtual
{
  public:

     //! Default constructor.
     IoVtk( std::string fname );

     //! Destructor
     virtual ~IoVtk();

     //! Reads the order parameter \f$\psi\f$ and the problem parameter list
     //! from a legacy VTK file into the arguments.
     virtual void read( std::vector<double_complex> *psi,
                        Teuchos::ParameterList      *problemParams );

     //! Writes the  order parameter \f$\psi\f$ and the problem parameter list
     //! into a legac VTK file.
     //! The data is written such that the lexicographical ordering of the
     //! nodes is preserved and the resulting file contains a state \f$\psi\f
     //! that can be viewed using standard tools.
     virtual void write( const std::vector<double_complex> &psi,
                         const Teuchos::ParameterList      &problemParams,
                         StaggeredGrid                     &sGrid          );

  protected:
  private:

      //! joins a vector of strings to one string with a separator string sep
      std::string strJoin( const std::vector<std::string> & vec,
                           const std::string              & sep  );

};