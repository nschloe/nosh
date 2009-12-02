/*
 * GridUniformVirtual.h
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schlšmer
 */

#ifndef GRIDUNIFORMVIRTUAL_H_
#define GRIDUNIFORMVIRTUAL_H_

#include "GridVirtual.h"

class GridUniformVirtual: virtual public GridVirtual
{
public:
    GridUniformVirtual( double scaling = 0.0,
                        double h = 0.0,
                        double gridDomainArea = 0.0,
                        int    numGridPoints = 0,
                        int    numBoundaryPoints = 0);

    virtual
    ~GridUniformVirtual();

    void
    setScaling( const double alpha );

    double
    getH() const; //!< Returns the uniform mesh size \f$h\f$.


  protected:
    double h_;

  private:
};

#endif /* GRIDUNIFORMVIRTUAL_H_ */
