#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/eigs.h>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  Eigen::SparseMatrix<double> L, M, A, Q, B;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  vector_area_matrix(F, A);
  Q = igl::repdiag(L, 2) - A;
  B = igl::repdiag(M, 2);
  
  // Use eigen-decomposition to find an appropriate minimizer
  Eigen::MatrixXd sU;
  Eigen::VectorXd sS;
  igl::eigs(Q, B, 3, igl::EigsType::EIGS_TYPE_SM, sU, sS);
  
  U = Eigen::MatrixXd(sU.col(2)); // the third eigenvector is the first non-trivial soln
  U.resize(V.rows(), 2);
  
  // Rotate U by a canonical rotation
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(U.transpose() * U, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = U * svd.matrixV(); // V = inv of V'
}
