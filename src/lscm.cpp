#include "lscm.h"
#include <igl/massmatrix.h>
#include <vector_area_matrix.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include <igl/repdiag.h>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int V_size = F.maxCoeff()+1;

  // Get vector area matrix
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F,A);

  // Get Laplacian
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  std::cout << L.row(0) << std::endl;

  // Create Q;
  Eigen::MatrixXd Q(V_size *2, V_size * 2);
  Q.block(0,0,V_size, V_size) = Eigen::MatrixXd(L);
  std::cout << "block one" << std::endl;
  Q.block(V_size,V_size, V_size, V_size) = L;
  std::cout << "block two" << std::endl;
  Eigen::SparseMatrix<double> Q_sparse = Q.sparseView() - A;
  std::cout << "block operations" << std::endl;

  // Get M - mass matrix
  std::cout << "creating M..." << std::endl;
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);

  // Create B
  std::cout << "creating B..." << std::endl;
  Eigen::MatrixXd B(2*V_size,2*V_size);
  B.setZero();
  B.block(0,0,V_size, V_size) = M;
  B.block(V_size,V_size, V_size, V_size) = M;
  Eigen::SparseMatrix<double> B_sparse = B.sparseView();

  // get second-smallest eigen-value
  Eigen::MatrixXd sU;
  Eigen::VectorXd sS;
  std::cout << "gets here" << std::endl;
  igl::eigs(Q_sparse,B_sparse,2,igl::EIGS_TYPE_SM,sU,sS);
  std::cout << sS << std::endl;
}
