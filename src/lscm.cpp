#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>
#include <igl/svd3x3.h>
#include <iostream>
#include <Eigen/SVD>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // compose Q
  // Q = (L 0; 0 L) - A
  Eigen::SparseMatrix<double> L, L_diag, Q, A;
  igl::cotmatrix(V, F, L);
  igl::repdiag(L, 2, L_diag);
  vector_area_matrix(F, A);
  Q = L_diag - A;

  // compose B
  // B = (M 0; 0 M)
  Eigen::SparseMatrix<double> M, B;
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  igl::repdiag(M, 2, B);

  // now with B and Q we can solve the generalized eigen value problem
  Eigen::MatrixXd eigenvecs;
  Eigen::VectorXd eigenvals;
  igl::eigs(Q, B, 3, igl::EigsType::EIGS_TYPE_SM, eigenvecs, eigenvals);
  Eigen::VectorXd nontrivial_vec = eigenvecs.col(2);

  // reformat the eigenvector (unvectorize) to form U
  int n = F.maxCoeff() + 1;
  U.resize(n, 2);
  U.col(0) = nontrivial_vec.head(n);
  U.col(1) = nontrivial_vec.tail(n);

  // do svd on U^T * U to get rotation
  Eigen::MatrixXd U_rot;
  U_rot = U.transpose() * U;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(U_rot, Eigen::ComputeFullU | Eigen::ComputeFullV);

  // apply rotation
  U = U * svd.matrixU();
}
