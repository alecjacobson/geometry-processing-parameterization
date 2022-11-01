#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

using namespace Eigen;
using namespace std;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int n = F.maxCoeff() + 1;
  A.resize(n * 2, n * 2);

  vector<vector<int>> b;
  igl::boundary_loop(F, b);

  // construct A
  // two points scenario
  //  0   0   0  0.5
  //  0   0 -0.5   0
  //  0 -0.5   0   0
  // 0.5   0   0   0
  vector<Triplet<double>> triplets;
  for (int i = 0; i < b.size(); i++) {
    vector<int> bd = b[i];
    for (int j = 0; j< bd.size(); j++) {
      int p1 = bd[j];
      int p2_index = (j + 1) % bd.size();
      int p2 = bd[p2_index];

      // ensure symmetry
      // A
      triplets.push_back(Triplet<double>(p1, p2 + n, 0.5));
      triplets.push_back(Triplet<double>(p2, p1 + n, -0.5));
      // A_T
      triplets.push_back(Triplet<double>(p2 + n, p1, 0.5));
      triplets.push_back(Triplet<double>(p1 + n, p2, -0.5));
    }
  }
  A.setFromTriplets(triplets.begin(), triplets.end());
  A = 0.5 * A;
}
