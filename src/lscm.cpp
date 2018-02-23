#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/eigs.h>
#include <igl/cotmatrix.h>
// Given a 3D mesh (`V`,`F`) with boundary compute a 2D parameterization that
// minimizes the "least squares conformal" energy:
//
// \\[
// ∫_Ω ‖ ∇v - (∇u)^⊥ ‖² dA,
// \\]
//
// where u and v are the unknown (output) coordinates in the parametric domain
// `U`.
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of triangle indices into V
// Outputs:
//   U  #U by 2 list of mesh UV parameterization coordinates
//
typedef Eigen::Triplet<double> tuple;
void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  U = V.leftCols(2);
  int nV = V.rows();

  Eigen::SparseMatrix<double> A(2*nV, 2*nV);
  vector_area_matrix(F, A);

  //Adapted from Adam Sturge's implementation
  //The L matrix constructed based on the formula in the README 
  //always causes igl::eigs() to crash
  //I can not figure out why TAT
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> L_diag(2*nV, 2*nV);
  igl::repdiag(L, 2, L_diag);
  Eigen::SparseMatrix<double> Q;
  Q = L_diag - A;

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  Eigen::SparseMatrix<double> B(2*nV, 2*nV);
  igl::repdiag(M, 2, B);

  //Adapted from Adam Sturge's implementation
  Eigen::VectorXd sU,sS;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, sU, sS);
  Eigen::VectorXd U_vstack = sU.col(2);
  U.col(0) = U_vstack.topRows(nV);
  U.col(1) = U_vstack.bottomRows(nV);
}