#include "lscm.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include "vector_area_matrix.h"
#include <igl/repdiag.h>
#include <Eigen/SVD>
#include <iostream>
void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code

  Eigen::SparseMatrix<double> A; //vector area matrix
  Eigen::SparseMatrix<double> M; //mass
  Eigen::SparseMatrix<double> B; //diag mass
  Eigen::SparseMatrix<double> L; //cot
  Eigen::SparseMatrix<double> Q; //diag cot - A
  vector_area_matrix(F, A);
  //std::cout << "here1" << std::endl;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::cotmatrix(V, F, L);

  igl::repdiag(M, 2, B);
  igl::repdiag(L, 2, Q);
  Q = Q - 2 * A;
  Eigen::MatrixXd sU;
  Eigen::VectorXd sS;
  //general eig solve
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, sU, sS);

  //create stack u v matrix
  int n = V.rows();
  Eigen::MatrixXd stackU(n, 2);
  stackU.col(0) = sU.block(0,2, n, 1); //get the 3rd smallest eig v
  stackU.col(1) = sU.block(n,2, n, 1);
  //std::cout << "here2" << std::endl;

  //rotate to matrixU base
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(stackU.transpose() * stackU, Eigen::ComputeThinU | Eigen::ComputeThinV);
  U = stackU * svd.matrixU();
  //std::cout << "here" << std::endl;
  //U = stackU; 
}
