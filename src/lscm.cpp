#include "lscm.h"
#include "vector_area_matrix.h"

#include <Eigen/SVD>
#include <igl/eigs.h>
#include <igl/repdiag.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    // Solve optimization as a generalized Eigen value problem
    //      min_U U'QU  subject to  U'BU = 1

    int n = V.rows(); 

    // Compute Q & B

    Eigen::SparseMatrix<double> A, L, Q;
    igl::cotmatrix(V, F, L);
    igl::repdiag(L, 2, Q);
    vector_area_matrix(F, A);
    Q = Q - A;

    Eigen::SparseMatrix<double> M, B;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    igl::repdiag(M, 2, B);

    // Solve

    Eigen::MatrixXd sU;
    Eigen::VectorXd sS;
    igl::eigs(Q, B, 4, igl::EIGS_TYPE_SM, sU, sS);

    // Somehow, first 2 eigen value ~e^{-13}
    //      3, 4, ... looks more reasonable
    U.resize(n, 2);
    U << sU.col(2).topRows(n), sU.col(2).bottomRows(n);

    // Find Canonical rotation V, using PCA with SVD on U'U
    //      Right singular vectors V are eigenvectors of U'U
    // https://en.wikipedia.org/wiki/Principal_component_analysis#Singular_value_decomposition

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose() * U,  Eigen::ComputeFullV);
    U = U * svd.matrixV();

}