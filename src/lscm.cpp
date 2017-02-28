#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <iostream>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Get cotangent Laplacian
  Eigen::SparseMatrix<double> L(V.rows(), V.rows());
  igl::cotmatrix(V, F, L);
  
  // Get mass matrix
  Eigen::SparseMatrix<double> M(V.rows(), V.rows());
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  
  // Get vector area matrix A
  Eigen::SparseMatrix<double> A(2*V.rows(), 2*V.rows());
  vector_area_matrix(F, A);
  
  // Calculate Q
  Eigen::SparseMatrix<double> Q(2*V.rows(), 2*V.rows());
  igl::repdiag(L, 2, Q);
  Q = -Q + A;
  
  // Calculate B
  Eigen::SparseMatrix<double> B(2*V.rows(), 2*V.rows());
  igl::repdiag(M, 2, B);
  
  // Solve for non-zero lambda s.t. U'BU=1 & QU = lambdaBU
  Eigen::MatrixXd sU(2*V.rows(), 3); 
  Eigen::VectorXd sS(3);
  eigs(Q, B, 3, igl::EigsType::EIGS_TYPE_SM, sU, sS);
  U.resize(V.rows(), 2);
  U.col(0) = sU.col(2).head(V.rows());
  U.col(1) = sU.col(2).tail(V.rows());
  
  // Canonical Rotation
  Eigen::MatrixXd cov = U.transpose() * U;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
  U = U * svd.matrixV().transpose();
}
