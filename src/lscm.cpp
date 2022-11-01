#include "lscm.h"
#include <igl/eigs.h>
#include <igl/repdiag.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include "vector_area_matrix.h"
#include <Eigen/SVD>

using namespace Eigen;
using namespace std;

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // construct matrix Q
  SparseMatrix<double> L, L_rep, Q, A, M, B;
  vector_area_matrix(F, A);
  igl::cotmatrix(V, F, L);
  igl::repdiag(L, 2, L_rep);
  Q = L_rep - A;

  // similarly, construct matrix B
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::repdiag(M, 2, B);

  // start to solve the generalized Eigen value problem
  MatrixXd eVector;
  VectorXd eValue;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eVector, eValue);
  // extract the Fiedler vector: 3rd -- first 2 trivial
  int n = V.rows();
  U.resize(n, 2);
  U.col(0) = eVector.col(2).block(0, 0, n, 1);
  U.col(1) = eVector.col(2).block(n, 0, n, 1);
  // using PCA: SVD solves
  JacobiSVD<MatrixXd> SVD(U.transpose() * U, ComputeThinU | ComputeThinV);
  U = U * SVD.matrixU();
}
