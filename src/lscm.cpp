#include "lscm.h"
#include "vector_area_matrix.h"
#include <Eigen/SVD>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/eigs.h>
#include <igl/repdiag.h>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = V.leftCols(2);
  Eigen::SparseMatrix<double> L,A,Q;
  igl::cotmatrix(V,F,L);
  vector_area_matrix(F,A);
  igl::repdiag(L,2,Q);
  Q = Q - A;
  Eigen::SparseMatrix<double> M,B;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::repdiag(M,2,B);

  Eigen::MatrixXd sU;
  Eigen::VectorXd sS;
  igl::eigs(Q,B,3,igl::EIGS_TYPE_SM,sU,sS);

  U.resize(V.rows(), 2);
  U.col(0) = sU.col(2).head(V.rows());
  U.col(1) = sU.col(2).tail(V.rows());
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose() * U, Eigen::ComputeThinU | Eigen::ComputeThinV);
  U = U * svd.matrixU();
  // since U is flipped due to boundary issue, I times -1 to flip the normal.
  U.col(0) *= -1;
}
