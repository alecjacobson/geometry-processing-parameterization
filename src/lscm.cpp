#include "lscm.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/repdiag.h"
#include "igl/eigs.h"
#include "Eigen/src/SVD/JacobiSVD.h"
#include "vector_area_matrix.h"
#include <iostream>

using namespace std;
using namespace igl;
using namespace Eigen;

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  U = V.leftCols(2);
  // Calculate L
  Eigen::SparseMatrix<double> L(V.rows(), V.rows());
  L.setZero();
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F, A);

  // Set up the quadratic matrices
  Eigen::SparseMatrix<double> L_rep;
  igl::repdiag(L, 2, L_rep);
  Eigen::SparseMatrix<double> Q = L_rep - A;
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  Eigen::SparseMatrix<double> B;
  igl::repdiag(M, 2, B);

  // Solving the generalized eigen vector problem
  Eigen::MatrixXd sU;
  Eigen::VectorXd sS;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, sU, sS);
  U.col(0) = sU.block(0, 2, V.rows(), 1);
  U.col(1) = sU.block(V.rows(), 2, V.rows(), 1);

  // Align the parameterization
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose() * U, Eigen::ComputeThinU | Eigen::ComputeThinV);
  U = U * svd.matrixV();
}
