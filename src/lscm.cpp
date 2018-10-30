#include "lscm.h"
#include <Eigen/SparseCore>
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>
#include <Eigen/SVD>
#include <Eigen/LU>

void lscm(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  int n = V.rows();
  U.resize(n, 2);

  Eigen::SparseMatrix<double> A, Q, L, L2;

  vector_area_matrix(F, A);
  igl::cotmatrix(V, F, L);
  igl::repdiag(L, 2, L2);
  Q = L2 - A;

  Eigen::SparseMatrix<double> M, B;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::repdiag(M, 2, B);

  // Solve general eigen value problem
  Eigen::MatrixXd eigenVec;
  Eigen::VectorXd eigenVal;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eigenVec, eigenVal);
  U.col(0) << eigenVec.col(2).head(n);
  U.col(1) << eigenVec.col(2).tail(n);

  // SVD
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose() * U, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = U * svd.matrixU();
}
