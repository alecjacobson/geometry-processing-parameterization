#include "lscm.h"
#include "igl/cotmatrix.h"
#include "vector_area_matrix.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>
#include <Eigen/SVD>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int n = V.rows();
  Eigen::SparseMatrix<double> L(n, n);
  Eigen::SparseMatrix<double> Q(n * 2, n * 2);
  Eigen::SparseMatrix<double> M(n, n);
  Eigen::SparseMatrix<double> M2(n * 2, n * 2);
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F, A);

  // std::vector<Eigen::Triplet<double> > tripletList;
  // for (int k = 0; k < L.outerSize(); ++k) {
  //   for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
  //     tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
  //     tripletList.push_back(Eigen::Triplet<double>(n + it.row(), n + it.col(), it.value()));
  //   }
  // }
  // Q.setFromTriplets(tripletList.begin(), tripletList.end());
  igl::repdiag(L, 2, Q);
  Q = Q - A;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
  igl::repdiag(M, 2, M2);
  Eigen::MatrixXd sU = Eigen::MatrixXd::Zero(n, 2);
  Eigen::VectorXd sS = Eigen::VectorXd::Zero(2);

  igl::eigs(Q, M2, 3, igl::EIGS_TYPE_SM, sU, sS);
  U = Eigen::MatrixXd::Zero(V.rows(), 2);
  U.col(0) = sU.block(0, 2, n, 1);
  U.col(1) = sU.block(n, 2, n, 1);

  Eigen::MatrixXd UT = U.transpose();
  Eigen::MatrixXd UTU = UT * U;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(UTU, Eigen::ComputeThinU | Eigen::ComputeThinV);
  U = U * svd.matrixU();
}
