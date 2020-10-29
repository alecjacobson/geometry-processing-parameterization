#include "lscm.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>
#include <igl/repdiag.h>
#include "vector_area_matrix.h"

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  U.resize(V.rows(), 2);
  // Quardratic Coefficients:
  Eigen::SparseMatrix<double> L, Q, A;
  igl::cotmatrix(V, F, L);
  vector_area_matrix(F, A);
  igl::repdiag(L, 2, Q);
  Q = Q - A;

  // Free Boundary:
  Eigen::SparseMatrix<double> M, B;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::repdiag(M, 2, B);
  Eigen::MatrixXd sU;
  Eigen::VectorXd sS;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, sU, sS);
  // First two are trivial solutions
  U.col(0) = sU.col(2).head(V.rows());
  U.col(1) = sU.col(2).tail(V.rows());

  // Apply SVD:
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose()*U, 
                    Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = U*svd.matrixU();
}
