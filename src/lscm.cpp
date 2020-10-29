#include "lscm.h"
#include "vector_area_matrix.h"
#include "igl/repdiag.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/eigs.h"
#include <Eigen/SVD>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{ 
  // Construct matrix Q = (L 0; 0 L) - A where L is the laplacian and A is the vector area
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> L_d;
  igl::repdiag(L, 2, L_d);
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F, A);
  Eigen::SparseMatrix<double> Q = L_d  - 2.0 * A;

  // Construct the matrix B = (M 0; 0 M) where M is the mass matrix
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_BARYCENTRIC, M);
  Eigen::SparseMatrix<double> B;
  igl::repdiag(M, 2, B);

  // Find parameterization coordinates by a eigenvector of Qv = lambda * B * v with v non-zero eigenvalue
  Eigen::VectorXd eigen_vals;
  Eigen::MatrixXd eigen_vecs;
  igl::eigs(Q, B, 5, igl::EigsType::EIGS_TYPE_SM, eigen_vecs, eigen_vals);
  
  int idx = 1; // the second eigenvector calculated should be the Fielder vector
   double tol = 1e-8; // Pick an eigenvector sufficiently away from 0
  // hope that an eigen_vec with non-zero eigen_val has been found among the calculated evals
  // otherwise just take the largest evec with the largest evec. 
  while (abs(eigen_vals[idx]) <= tol) {
	idx++;
  }
  if (idx >= 5) {
	idx = 4;
  }

  // Reshape eigenvector into U
  int n_v = F.maxCoeff() + 1;
  U.resize(n_v, 2);
  U.col(0) = eigen_vecs.col(idx).head(n_v);
  U.col(1) = eigen_vecs.col(idx).tail(n_v);
  
  // Do SVD on coordinates U to find a rotation on the coordinates U
  Eigen::Matrix2d covar = U.transpose() * U;
  Eigen::JacobiSVD<Eigen::Matrix2d> svd(covar, Eigen::ComputeFullU);
  for (int i = 0; i < U.rows(); i++) {
	U.row(i) = U.row(i) * svd.matrixU(); 
  }
}
