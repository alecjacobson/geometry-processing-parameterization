#include "lscm.h"
#include "vector_area_matrix.h"

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>

#include <assert.h>


void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    const int n = V.rows();
    const int n2 = F.maxCoeff() + 1;

    assert( n == n2 );
    
    Eigen::SparseMatrix< double > M( n, n );
    igl::massmatrix( V, F, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI, M );

    // make matrix B such that
    //
    //       M 0
    //  B =  0 M
    //
    Eigen::SparseMatrix< double > B( 2*n, 2*n );
    igl::repdiag( M, 2, B );

    Eigen::SparseMatrix< double > L( n, n );
    igl::cotmatrix( V, F, L );

    Eigen::SparseMatrix< double > A( 2*n, 2*n );
    vector_area_matrix( F, A );

    
    Eigen::SparseMatrix< double > Q( 2*n, 2*n );
    Eigen::SparseMatrix< double > LL = (-L)*.5;
    igl::repdiag( LL, 2, Q );
    Q = Q - A;

    // throw into eigs, choose 3rd eigenvector, Fiedler
    // not to be confused with the Arthur Fielder Footbridge in Boston.
    const int k = 3;
    Eigen::MatrixXd sU( 2*n, k );
    Eigen::VectorXd sS( k );
    igl::eigs( Q, B, 3, igl::EigsType::EIGS_TYPE_SM, sU, sS );

    U.resize( n, 2 );
    U.col( 0 ) = sU.col(2).block( 0, 0, n, 1 );
    U.col( 1 ) = sU.col(2).block( n, 0, n, 1 );

    // from our closest_rotation in registration....
    Eigen::JacobiSVD< Eigen::MatrixXd > svd( U, Eigen::ComputeThinU | Eigen::ComputeThinV );

    U = U * svd.matrixV();
}
