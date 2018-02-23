#include "lscm.h"
#include "vector_area_matrix.h"

#include "igl/eigs.h"
#include "igl/repdiag.h"
#include "igl/massmatrix.h"

#include <iostream>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	Eigen::SparseMatrix<double> L, A, M;
  igl::cotmatrix(V, F, L); 
  vector_area_matrix(F, A);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  Eigen::SparseMatrix<double> Q, B;
  igl::repdiag(M, 2, B);
  igl::repdiag(L, 2, Q);
  Q -= A;

  Eigen::MatrixXd eig_vecs;
  Eigen::VectorXd eig_vals;
  igl::eigs(Q, B, 5, igl::EIGS_TYPE_SM, eig_vecs, eig_vals);

  U = Eigen::Map<Eigen::MatrixXd>(eig_vecs.col(1).data(), V.rows(), 2);

  // Canonical rotation
  Eigen::MatrixXd UTU = U.transpose() * U;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(UTU, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = U * svd.matrixU().transpose();
}
