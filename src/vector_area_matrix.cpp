#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  typedef Eigen::Triplet<double> T;
  
  int V_size = F.maxCoeff()+1;
  A.resize(V_size*2,V_size*2);

  // Construct A_bar
  std::vector<T> tripletList;
  std::vector<std::vector<int>> bloop;
  igl::boundary_loop(F, bloop);
  for (int a = 0; a < bloop.size(); ++a) {
    for (int b = 0; b < bloop[a].size(); ++b) {
      const int& i = bloop[a][b];
      const int& j = bloop[a][(b + 1)%bloop[a].size()];
      
      int ix = i;
      int iy = V_size + i;
      int jx = j;
      int jy = V_size + j;

      // determinant ix*jy - iy*jx
      // ix*jy: row -> position of jy, col -> position of ix
      tripletList.push_back(T(jy, ix, 0.5));

      // iy*jx: row -> position of jx, col -> position of iy
      tripletList.push_back(T(jx, iy, -0.5));
    } // end loop b
  } // end loop a

  Eigen::SparseMatrix<double> A_bar(V_size*2, V_size*2);
  A_bar.setFromTriplets(tripletList.begin(), tripletList.end());

  A = 0.5*(A_bar + Eigen::SparseMatrix<double>(A_bar.transpose()));
}

