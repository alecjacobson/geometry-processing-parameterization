#include "vector_area_matrix.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int n = F.maxCoeff()+1;
  A.resize(n*2,n*2);

  std::vector<std::vector<int>> bdry;
  igl::boundary_loop(F,bdry);
  // Outputs:
  // L  list of loops where L[i] = ordered list of boundary vertices in loop i

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  for (int i=0; i < bdry.size(); i++){
    for (int j=0; j < bdry[i].size(); j++){
      int pi = bdry[i][j];
      int pj = bdry[i][(j+1)%bdry[i].size()];
    
      tripletList.push_back(T(pi, pj + n, 0.5));
      tripletList.push_back(T(pi + n, pj, -0.5));
    }
  }

  A.setFromTriplets(tripletList.begin(), tripletList.end());
  Eigen::SparseMatrix<double> AT;
  AT.resize(n*2, n*2);
  AT = A.transpose();
  A = 0.5 * (A + AT); // symmetrize
}
