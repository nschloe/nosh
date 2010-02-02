/*
 * GridUniformVirtual.h
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schl\omer
 */

#ifndef GRIDUNIFORMVIRTUAL_H_
#define GRIDUNIFORMVIRTUAL_H_

#include "GridVirtual.h"

class GridUniformVirtual: virtual public GridVirtual
{
public:
    //! Default constructor.
    GridUniformVirtual( double scaling,
                        double h,
                        double gridDomainArea,
                        unsigned int numGridPoints,
                        unsigned int numBoundaryPoints);

    //! Class constructor that only initializes the data members of this class.
    GridUniformVirtual();

    virtual
    ~GridUniformVirtual();

    double
    getUniformH() const; //!< Returns the uniform mesh size \f$h\f$.

  protected:
  private:
};

#endif /* GRIDUNIFORMVIRTUAL_H_ */
