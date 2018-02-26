#include "lscm.h"
#include <igl/repdiag.h>
#include <vector_area_matrix.h>
#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <Eigen/SVD>
#include <igl/massmatrix.h>



void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> M_diag;
  Eigen::SparseMatrix<double> L_diag;
  Eigen::MatrixXd U_temp;
  Eigen::VectorXd S;
  
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  vector_area_matrix(F, A);
  igl::repdiag(L, 2 , L_diag);
  igl::repdiag(M, 2 , M_diag);
  A = L_diag - A;
  igl::eigs(A, M_diag, 5, igl::EIGS_TYPE_SM, U_temp, S);
  
  U.resize(V.rows(), 2);
  U.col(0) = U_temp.col(1).head(V.rows());
  U.col(1) = U_temp.col(1).tail(V.rows());
  Eigen::JacobiSVD<Eigen::Matrix2d> svd(U.transpose() * U, Eigen::ComputeFullU | Eigen::ComputeFullV);

  U *= svd.matrixV();
}