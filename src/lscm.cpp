#include "lscm.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include "vector_area_matrix.h"
#include <igl/repdiag.h>
#include <igl/eigs.h>
#include <iostream>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  U = V.leftCols(2);

  std::cout << "LSCM start" << std::endl;
  int dim = V.rows();


  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);

  Eigen::SparseMatrix<double> VA;
  vector_area_matrix(F,VA);

  Eigen::SparseMatrix<double> LD(2*dim,2*dim);
  igl::repdiag(L,2,LD);

  Eigen::SparseMatrix<double> A = LD + VA;

  Eigen::SparseMatrix<double> B(2*dim,2*dim);
  igl::repdiag(M,2,B);

  Eigen::MatrixXd EVec;
  Eigen::VectorXd EVal;

  /*
  std::cout <<A << std::endl;
  std::cout <<B << std::endl;
  */

  std::cout << "Eigs?" << std::endl;
  igl::eigs(A,B,3,igl::EIGS_TYPE_SM,EVec,EVal);
  std::cout << EVal.transpose() << std::endl;
  std::cout << "LSCM end" << std::endl;
  U.col(0) = EVec.col(2).head(dim);
  U.col(1) = EVec.col(2).tail(dim);

}
