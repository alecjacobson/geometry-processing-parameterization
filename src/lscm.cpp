#include "lscm.h"

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int n = F.maxCoeff() + 1;

  Eigen::SparseMatrix<double> A, L, Q, M, B;

  vector_area_matrix(F, A); // area matrix
  igl::cotmatrix(V, F, L); // cot matrix
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M); // mass matrix

  igl::repdiag(M, 2, B);
  igl::repdiag(L, 2, Q);
  Q = Q - 2*A;

  Eigen::MatrixXd eigvec;
  Eigen::VectorXd eigval;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eigvec, eigval);
  // after 2 trivial solutions

  U.resize(n, 2);
  U.col(0) = eigvec.col(2).head(n); // u
  U.col(1) = eigvec.col(2).tail(n); // v

  // canonical rotation
  Eigen::MatrixXd tmp = U.transpose() * U;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(tmp, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = U * svd.matrixU();
}