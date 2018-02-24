#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  A.resize(2*F.maxCoeff()+2,2*F.maxCoeff()+2);
  Eigen::VectorXi B;
  igl::boundary_loop(F, B);
  int n = F.maxCoeff() + 1;
  for (int i = 0;i<B.rows() - 1;i++) {
  	A.coeffRef(B(i), n + B(i+1)) = 0.25;
  	A.coeffRef(B(i+1), n + B(i)) = -0.25;
  }
  Eigen::SparseMatrix<double> AT = A.transpose();
  A = A + AT;
}

