/*
 * GridUniformVirtual.h
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef GRID_UNIFORM_ABSTRACT_H_
#define GRID_UNIFORM_ABSTRACT_H_

#include "Recti_Grid_Abstract.h"

namespace Recti
{
  namespace Grid
  {

class UniformAbstract:
            virtual public Abstract
{
public:
    //! Default constructor.
    UniformAbstract ( double h,
                      double gridDomainArea,
                      unsigned int numGridPoints );

    //! Class constructor that only initializes the data members of this class.
    UniformAbstract();

    virtual
    ~UniformAbstract();

    double
    getUniformH() const; //!< Returns the uniform mesh size \f$h\f$.

protected:
private:
};

  } // namespace Grid
} // namespace Recti

#endif /* GRID_UNIFORM_ABSTRACT_H_ */
