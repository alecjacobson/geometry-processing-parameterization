#include "vector_area_matrix.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  A.resize(F.maxCoeff()+1,F.maxCoeff()+1);
}

