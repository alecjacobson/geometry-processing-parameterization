#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/eigs.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <Eigen/SVD>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U) {
    
    // construct the cotangent Laplacian
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);
    
    // construct the vector area matrix
    Eigen::SparseMatrix<double> A;
    vector_area_matrix(F, A);
    
    // construct the Q matrix in our eigenvalue problem using a matrix with L repeated
    // along the diagonal and our symmetrized A.
    Eigen::SparseMatrix<double> Q;
    igl::repdiag(L, 2, Q);
    Q -= A;
    
    // construct the mass matrix
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    
    // construct the B matrix in our eigenvalue problem by repeating M twice along
    // the diagonal
    Eigen::SparseMatrix<double> B;
    igl::repdiag(M, 2, B);
    
    // compute 3 solution eigenpairs using EIGS_TYPE_SM (we need 3 because the first 2
    // are trivial solutions)
    Eigen::MatrixXd sU;
    Eigen::VectorXd sS;
    igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, sU, sS);
    
    // select the smallest non-zero eigenvector and place the corresponding u coordinates
    // in the first column of U and the v coordinates in the second column
    Eigen::VectorXd eig = sU.col(2);
    int n = V.rows();
    U.resize(n, 2);
    U.col(0) = eig.topRows(n);
    U.col(1) = eig.bottomRows(n);
    
    // compute the singular value decomposition of the covariance matrix to obtain
    // the canonical rotation
    Eigen::MatrixXd cov = U.transpose()*U;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(cov, Eigen::ComputeFullV);
    Eigen::MatrixXd rot = svd.matrixV();
    
    // rotate U by the canonical rotation
    U *= rot;   
}
