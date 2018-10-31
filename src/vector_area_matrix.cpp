#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
#include <vector>
#include <iostream>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // U^T * A * U
  // compose A to apply 1/2 Σ |u_i u_j| over all boundary edges
  // given that U is in a vectorized form

  // first grab all the boundary loops
  std::vector<std::vector<int>> boundaries;
  igl::boundary_loop(F, boundaries);

  // now loop over them to compose A
  int n = F.maxCoeff() + 1;
  int u_i, u_j;
  std::vector<Eigen::Triplet<double>> triplets;
  for (const auto& boundary: boundaries) {
    for (int i = 0; i < boundary.size(); i++) {
      u_i = boundary[i];
      u_j = boundary[(i + 1) % boundary.size()];
      // |u_i u_j| = u_i(0) * u_j(1) - u_i(1) * u_j(0)
      // effectively, u_x(1) = u_x(0) + n, since U is vectorized
      triplets.push_back(Eigen::Triplet<double>(u_i, u_j + n, 0.25));
      triplets.push_back(Eigen::Triplet<double>(u_i + n, u_j, -0.25));
    }
  }
  Eigen::SparseMatrix<double> A_asym(2 * n, 2 * n);
  A_asym.setFromTriplets(triplets.begin(), triplets.end());

  // we want a symmetric matrix -> (1/2) (A + A') 
  // (applied the 1/2 earlier)
  Eigen::SparseMatrix<double> A_asym_T(2 * n, 2 * n);
  A_asym_T = A_asym.transpose();
  A.resize(2 * n, 2 * n);
  A = A_asym + A_asym_T;
}

