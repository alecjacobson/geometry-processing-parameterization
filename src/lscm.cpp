#include "lscm.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include "vector_area_matrix.h"
#include <igl/repdiag.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/min_quad_with_fixed.h>
#include <igl/boundary_loop.h>
#include <igl/eigs.h>
#include <cmath>
#include <igl/polar_svd.h>
#include <iostream>


void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code

  // get the mass matrix
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

  // get the cotmatrix
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // vectorized version of U
  Eigen::MatrixXd U_vec;

  // get the area vector A
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F, A);

  // get the diag repeated matrix of L
  Eigen::SparseMatrix<double> L_rep;
  igl::repdiag(L, 2, L_rep);

  // get the quadratic coefficients out
  Eigen::SparseMatrix<double> Q = L_rep - A;


  // get the diag repeated matrix of M
  Eigen::SparseMatrix<double> M_rep;
  igl::repdiag(M, 2, M_rep);

  /*
  Need to minimize 0.5*U'*Q*U wrt U such that U'*M_rep*U=1.

  Solve the Generalized Eigenvalue problem: Q*U = lambda*M_rep*U.
  But we also have U'*M_rep*U=1.

  Therefore, putting that in prev eqn, we have lambda=U'*Q*U.
  So, minimizing the original problem is nothing but chosing the eigenvector with the 
  least eigenvalue!
  */

  Eigen::VectorXd lambda;

  // get the smallest eigen values
  // but the smallest one seems to be very close to zero
  // take the smallest but nontrivial (nonzero) eigen value
  igl::eigs(Q, M_rep, 5, igl::EIGS_TYPE_SM, U_vec, lambda);
 

  double best_lambda = lambda(lambda.size()-1);
  int best_lamda_idx = 5;
  for (int i = 0; i<lambda.size(); i++) {
    if (std::abs(lambda(i)) > 0.00000001) {
      best_lambda = lambda(i);
      best_lamda_idx = i;
      break; // pick the first non trivial eigen and quit loop
    }
  }
  
  // print the eigen values to check
  // std::cout<<lambda<<std::endl;
  //  4.76244e-13 
  //  3.43539e-13
  // -0.000837413
  // -0.000837413
  //   -0.0149227
  //  4.76244e-13
  //  3.43539e-13
  // -0.000837413
  // -0.000837413
  // -0.0149227

  // std::cout<<"Eigen value used: "<<best_lamda_idx<<" th."<<std::endl;

  // pick the smallest but non zero eigen value

  Eigen::MatrixXd U_vec_best(2*V.rows(), 1);
  U_vec_best = U_vec.col(best_lamda_idx);

  double lambda_out = (U_vec_best.transpose()*Q*U_vec_best*1.0)(0, 0);


  // rescale such that U_vec'*Q*U_vec = lambda
  U_vec_best /= std::sqrt(std::abs(best_lambda/lambda_out));


  // reshape U_vec to n X 2
  Eigen::MatrixXd U_bef_rota(V.rows(), 2);
  U_bef_rota.block(0, 0, V.rows(), 1) = U_vec_best.block(0, 0, V.rows(), 1);
  U_bef_rota.block(0, 1, V.rows(), 1) = U_vec_best.block(V.rows(), 0, V.rows(), 1);


  // rotation and translation
  Eigen::MatrixXd R(2,2);
  Eigen::MatrixXd T(2,2);


  Eigen::MatrixXd U_corr = U_bef_rota.transpose()*U_bef_rota;

  // get the canonical rotation 
  igl::polar_svd(U_corr, R, T);

  U.resize(V.rows(), 2);


  // apply rotation
  for (int row = 0; row<U_bef_rota.rows(); row++) {
    U.row(row) = (R * U_bef_rota.row(row).transpose()).transpose();
  }
  
}
