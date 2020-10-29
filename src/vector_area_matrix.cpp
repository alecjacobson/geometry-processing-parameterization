#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int n = F.maxCoeff() + 1;
  
  std::vector<std::vector<int>> bds;
  igl::boundary_loop(F, bds);
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  for (auto bd : bds) {
    int bd_len = bd.size();
    for (int vi = 0; vi < bd_len; vi++) {
      int i = bd[vi];
      int j = bd[(vi+1)%bd_len];
      tripletList.push_back(T(i, n + j, 1));
      tripletList.push_back(T(n + i, j, -1));
      // ensure that A is symmetric (add values for A^T)
      tripletList.push_back(T(n + j, i, 1));
      tripletList.push_back(T(j, n + i, -1));
    }
  }
  A.resize(2*n, 2*n);
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  A *= .5; // I'm not sure why 1/2 is the coeff. that worked; should be (1/2)*(1/2)?
}

