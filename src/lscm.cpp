#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/eigs.h>
#include <igl/repdiag.h>
#include <Eigen/src/SVD/JacobiSVD.h>
#include <vector>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int n = V.rows();
  Eigen::SparseMatrix<double> laplacian(n,n);
  igl::cotmatrix(V, F, laplacian);

  Eigen::SparseMatrix<double> mass(n,n);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

  Eigen::SparseMatrix<double> area;
  vector_area_matrix(F, area);

  Eigen::SparseMatrix<double> A = area;
  Eigen::SparseMatrix<double> Q(2 * n, 2 * n);
  igl::repdiag(laplacian, 2, Q);
  Q = Q - A;

  Eigen::SparseMatrix<double> B(2 * n, 2 * n);
  igl::repdiag(mass, 2, B);

  Eigen::MatrixXd eig_vecs;
  Eigen::VectorXd eig_vals;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eig_vecs, eig_vals);

  U.resize(n, 2);
  U.col(0) = eig_vecs.block(0, 2, n, 1);
  U.col(1) = eig_vecs.block(n, 2, n, 1);

  Eigen::Matrix2d k = U.transpose()*U;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(k, Eigen::ComputeThinU | Eigen::ComputeThinV);
  U = U*svd.matrixU();

}
