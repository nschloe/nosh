// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2010, 2011  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER

#ifndef GINLA_KEOFACTORY_H
#define GINLA_KEOFACTORY_H
// =============================================================================
// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include "Ginla_config.h"

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif
#include <Teuchos_Tuple.hpp>

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#include "Ginla_MagneticVectorPotential.hpp"
// =============================================================================
// forward declarations
namespace Ginla {
    class StkMesh;
}
// =============================================================================
namespace Ginla {
// =============================================================================
class KeoFactory
{
public:
    enum EMatrixType { MATRIX_TYPE_REGULAR,
                       MATRIX_TYPE_DMU
                     };

public:
    KeoFactory( const Teuchos::RCP<const Ginla::StkMesh>      & mesh,
                const Teuchos::RCP<const Epetra_Vector>            & thickness,
                const Teuchos::RCP<Ginla::MagneticVectorPotential> & mvp
              );

    // Destructor.
    ~KeoFactory();

    const Teuchos::RCP<const Ginla::StkMesh>
    getMesh() const;

    void
    updateParameters( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams
                    ) const;

    const Teuchos::RCP<const LOCA::ParameterVector>
    getMvpParameters() const;

    Teuchos::RCP<const Epetra_CrsGraph>
    getKeoGraph() const;

    Teuchos::RCP<const Epetra_CrsMatrix>
    getKeo() const;

    Teuchos::RCP<const Epetra_CrsMatrix>
    getKeoDMu() const;

protected:

private:
    const Teuchos::RCP<Epetra_CrsGraph>
    buildKeoGraph_() const;

    Teuchos::RCP<Epetra_CrsMatrix>
    buildKeo_( const EMatrixType matrixType ) const;

    void
    fillKeo_( const Teuchos::RCP<Epetra_CrsMatrix> & keoMatrix,
              const EMatrixType matrixType
            ) const;

    void
    fillKeoCellEdges_( const Teuchos::RCP<Epetra_CrsMatrix> & keoMatrix,
                       const EMatrixType matrixType,
                       const Teuchos::ArrayRCP<const DoubleVector> & edgeCoefficientsFallback
                     ) const;

    void
    fillKeoEdges_( const Teuchos::RCP<Epetra_CrsMatrix> & keoMatrix,
                   const EMatrixType matrixType,
                   const Teuchos::ArrayRCP<const double> & edgeCoefficients
                 ) const;

private:
#ifdef GINLA_TEUCHOS_TIME_MONITOR
    const Teuchos::RCP<Teuchos::Time> buildKeoTime_;
    const Teuchos::RCP<Teuchos::Time> buildKeoGraphTime_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime1_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime2_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime3_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime4_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime5_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime6_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime7_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime8_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime9_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime10_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime11_;
    const Teuchos::RCP<Teuchos::Time> buildKeoTime12_;
#endif
    const Teuchos::RCP<const Ginla::StkMesh> mesh_;
    const Teuchos::RCP<const Epetra_Vector> thickness_;
    const Teuchos::RCP<Ginla::MagneticVectorPotential> mvp_;
    const Teuchos::RCP<const Epetra_CrsGraph> keoGraph_;
    const Teuchos::RCP<Epetra_CrsMatrix> keo_;
    mutable Teuchos::RCP<LOCA::ParameterVector> keoBuildParameters_;
    const Teuchos::RCP<Epetra_CrsMatrix> keoDMu_;
    mutable Teuchos::RCP<LOCA::ParameterVector> keoDMuBuildParameters_;
};
// =============================================================================
} // namespace Ginla

#endif // GINLA_KEOFACTORY_H
