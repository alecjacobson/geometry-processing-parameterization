#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/eigs.h>
#include <Eigen/Eigenvalues>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  U = V.leftCols(2);

  // Compute cotmatrix
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  L *= -1.0;  // positive semi-definite

  // Compute vector area matrix
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F, A);

  // Construct Q
  Eigen::SparseMatrix<double> Q = 0.5*igl::repdiag(L, 2) - A;

  // Construct B
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
  Eigen::SparseMatrix<double> B = igl::repdiag(M, 2);

  // Solve
  Eigen::MatrixXd eig_vectors;
  Eigen::VectorXd eig_values;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eig_vectors, eig_values);

  // Convert to (u,v)
  U.block(0, 0, V.rows(), 1) = eig_vectors.block(0, 0, V.rows(), 1);
  U.block(0, 1, V.rows(), 1) = eig_vectors.block(V.rows(), 0, V.rows(), 1);

  // PCA
  Eigen::MatrixXd cov = U.transpose()*U;
  Eigen::EigenSolver<Eigen::MatrixXd> solver(cov);
  Eigen::MatrixXd rot = solver.eigenvectors().real();

  U = U*rot;
}
