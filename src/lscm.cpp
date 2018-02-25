#include "lscm.h"
#include <igl/massmatrix.h>
#include <vector_area_matrix.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include <igl/repdiag.h>
#include <Eigen/SVD>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int V_size = F.maxCoeff()+1;

  // Get vector area matrix
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F,A);

  // Get Laplacian
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  // Create Q;
  Eigen::SparseMatrix<double> Q;
  igl::repdiag(L, 2, Q);
  Q = Q - A;

  // Get M - mass matrix
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  // Create B
  Eigen::SparseMatrix<double> B;
  igl::repdiag(M, 2, B);

  // get second-smallest eigen-value
  // minv½v𝖳Qv subject to v𝖳Bv=1
  // v is 1D with 2*V_size entries, eigenvectors also have this size
  Eigen::MatrixXd sU;
  Eigen::VectorXd sS;
  int k = 4; // eigs works in decreasing order, k = 4 works well for beetle
  igl::eigs(Q,B,k,igl::EIGS_TYPE_SM,sU,sS);

  // Set U from the eigenvector corresponding to smallest eigenvalue
  U.resize(V_size, 2);
  U.col(0) = sU.col(k-1).head(V_size); // x coordinates in first half
  U.col(1) = sU.col(k-1).tail(V_size); // y coordinates in second half

  //Canonical Rotation
  Eigen::MatrixXd covariance(2,2);
  covariance = U.transpose() * U;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd( covariance, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // matrixU and matrixV.transpose() are the same..
  // pick one and rotate U by it to align with x-axis
  U *= svd.matrixU();
}
