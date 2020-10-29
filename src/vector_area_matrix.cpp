#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
    const Eigen::MatrixXi & F,
    Eigen::SparseMatrix<double>& A)
{
  int n = F.maxCoeff()+1;
  A.resize(n * 2,n * 2);

  std::vector<std::vector<int>> bL;
  igl::boundary_loop(F, bL);

  for (int x = 0; x < bL.size(); x++) {
    std::vector<int> b = bL[x];
    for (int y = 0; y < b.size(); y++) {
      int i = b[y];
      int j = b[(y + 1) % b.size()];
      A.coeffRef(i, n + j) = 0.5;
      A.coeffRef(j, n + i) = -0.5;
    }
  }

  // A to be symmetric
  Eigen::SparseMatrix<double> AT = A.transpose();
  A = A + AT;
  // However A *= 0.5 gives strange pattern.
}

