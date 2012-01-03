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

#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>

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
                       MATRIX_TYPE_DMU,
                       MATRIX_TYPE_DTHETA
                     };

public:
    KeoFactory( const Teuchos::RCP<const Ginla::StkMesh>           & mesh,
                const Teuchos::RCP<const Epetra_Vector>            & thickness,
                const Teuchos::RCP<Ginla::MagneticVectorPotential> & mvp
              );

    // Destructor.
    ~KeoFactory();

    const Epetra_Comm &
    getComm() const;

    void
    updateParameters( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams
                    ) const;

    const Teuchos::RCP<const LOCA::ParameterVector>
    getMvpParameters() const;

    Teuchos::RCP<const Epetra_FECrsMatrix>
    getKeo() const;

    Teuchos::RCP<const Epetra_FECrsMatrix>
    getKeoDMu() const;

    Teuchos::RCP<const Epetra_FECrsMatrix>
    getKeoDTheta() const;

protected:

private:
    const Teuchos::RCP<Epetra_FECrsGraph>
    buildKeoGraph_() const;

    void
    fillKeo_( const Teuchos::RCP<Epetra_FECrsMatrix> & keoMatrix,
              const EMatrixType matrixType
            ) const;

    void
    buildKeoGraphEdges_( const Teuchos::RCP<Epetra_FECrsGraph> & keoGraph ) const;

    void
    fillKeoEdges_( const Teuchos::RCP<Epetra_FECrsMatrix> & keoMatrix,
                   const EMatrixType matrixType,
                   const Teuchos::ArrayRCP<const double> & edgeCoefficients
                 ) const;

    void
    buildGlobalIndexCache_( const std::vector<stk::mesh::Entity*> & edges ) const;

    void
    buildAlphaCache_( const std::vector<stk::mesh::Entity*> & edges,
                      const Teuchos::ArrayRCP<const double> & edgeCoefficients
                    ) const;

    void
    buildKeoGraphCellEdges_( const Teuchos::RCP<Epetra_FECrsGraph> & keoGraph ) const;

    void
    fillKeoCellEdges_( const Teuchos::RCP<Epetra_FECrsMatrix> & keoMatrix,
                       const EMatrixType matrixType,
                       const Teuchos::ArrayRCP<const DoubleVector> & edgeCoefficientsFallback
                     ) const;

    void
    buildGlobalIndexFallbackCache_( const std::vector<stk::mesh::Entity*> & cells ) const;

    void
    buildAlphaFallbackCache_( const std::vector<stk::mesh::Entity*> & cells,
                              const Teuchos::ArrayRCP<const DoubleVector> & edgeCoefficientsFallback
                            ) const;

private:
#ifdef GINLA_TEUCHOS_TIME_MONITOR
    const Teuchos::RCP<Teuchos::Time> keoFillTime_;
    const Teuchos::RCP<Teuchos::Time> buildKeoGraphTime_;
#endif
    const Teuchos::RCP<const Ginla::StkMesh> mesh_;
    const Teuchos::RCP<const Epetra_Vector> thickness_;
    const Teuchos::RCP<Ginla::MagneticVectorPotential> mvp_;

    mutable Teuchos::ArrayRCP<Epetra_IntSerialDenseVector> globalIndexCache_;
    mutable bool globalIndexCacheUpToDate_;
    mutable Teuchos::ArrayRCP<Teuchos::ArrayRCP<Epetra_IntSerialDenseVector> > globalIndexFallbackCache_;
    mutable bool globalIndexFallbackCacheUpToDate_;

    const Teuchos::RCP<const Epetra_FECrsGraph> keoGraph_;
    const Teuchos::RCP<Epetra_FECrsMatrix> keo_;
    mutable Teuchos::RCP<LOCA::ParameterVector> keoBuildParameters_;
    const Teuchos::RCP<Epetra_FECrsMatrix> keoDMu_;
    mutable Teuchos::RCP<LOCA::ParameterVector> keoDMuBuildParameters_;
    const Teuchos::RCP<Epetra_FECrsMatrix> keoDTheta_;
    mutable Teuchos::RCP<LOCA::ParameterVector> keoDThetaBuildParameters_;

    mutable Teuchos::ArrayRCP<double> alphaCache_;
    mutable bool alphaCacheUpToDate_;
    mutable Teuchos::ArrayRCP<DoubleVector> alphaFallbackCache_;
    mutable bool alphaFallbackCacheUpToDate_;
};
// =============================================================================
} // namespace Ginla

#endif // GINLA_KEOFACTORY_H
