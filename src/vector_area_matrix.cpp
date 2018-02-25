#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
#include <vector>
#include <iostream>
void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  //A.resize(F.maxCoeff()+1,F.maxCoeff()+1);
  int n = F.maxCoeff()+1;
  A.resize(2 * n, 2 * n);

  std::vector<std::vector<int>> boundary_loop;
  igl::boundary_loop(F, boundary_loop);
  std::cout << boundary_loop.size() << std::endl;
  //follow the all boundary loop edge
  for (int b = 0; b < boundary_loop.size(); ++b){
    std::vector<int> Lb = boundary_loop[b];
    for (int i = 0; i < Lb.size(); ++i){
        int u_i = Lb[i];
        int u_j = Lb[(i+1) % Lb.size()];
        //u_j + n = v_j
        A.coeffRef(u_i, u_j + n) = 1;
        A.coeffRef(u_i+n, u_j) = -1;
    }
  }

  A = 0.5 * A;
  //symmetric A
  A = 0.5 * A + 0.5 * Eigen::SparseMatrix<double>(A.transpose());
}

