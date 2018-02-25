#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int V = F.maxCoeff()+1;
  A.resize(2*V,2*V);

  // get ordered vertices of boundary
  Eigen::VectorXi bnd;
  igl::boundary_loop(F,bnd);

  // loop over boundary vertices
  std::vector< Eigen::Triplet<double> > tripletList;
  int bnd_size = bnd.rows();
  for(int bnd_idx = 0; bnd_idx < bnd_size; bnd_idx += 1){
    int i_idx = bnd[bnd_idx];
    // get i's neighbor
    int j_idx = bnd[(bnd_idx + 1 )% bnd_size];
    // (ui_x * uj_y)  - (u_jx * ui_y)
    tripletList.push_back(Eigen::Triplet<double>(i_idx,j_idx + V, 0.5));
    tripletList.push_back(Eigen::Triplet<double>(i_idx + V,j_idx, -0.5));
  }
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  // Symmetrize A - take average of unsymmetrize A
  A = Eigen::SparseMatrix<double>(A.transpose()) + A;
  A *= 0.5;
}
