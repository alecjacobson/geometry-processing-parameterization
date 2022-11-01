#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
#include <vector>


void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  int V_size = F.maxCoeff()+1;
  A.resize(V_size*2,V_size*2);

  /*
  We need A to be such that U'AU gives the sum of signed areas where U is of dim 2*V_size X 1, i.e.,
  U(0:V_size) corresponds to first dimension and U(V_size:2*V_size) corresponds to the second dimension of the 2D parametrization.

  So, the u- and v- coordinate for the parametrized vertex indexed by i can be given by
  U(i) and U(i + V_size) respectively in the 2*V_size X 1 matrix notation.

  Let W be the matrix U represented as V_size X 2 (usual notation).
  Then the signed area for the edges i, j corresponds to the determinant:

  W(i, 0) W(j, 1) - W(i, 1) W(j, 0) which is same as 
  U(i) U(j + V_size) - U(i + V_size) U(j) in the 2*V_size X 1 notation.

  So, the matrix definition of A should be as follows:

  For all edges (i, j) on each boundary loop.
    A[i, j + V_size] = 1,
    A[i + V_size, j] = -1.

  Then symmetrize A. 
  That will be same as:
    A[i, j + V_size] = 1/2,
    A[i + V_size, j] = -1/2,
    A[j, i + V_size] = -1/2,
    A[j + V_size, i] = 1/2.
  */

 std::vector<Eigen::Triplet<double> > triplets;

 std::vector<std::vector<int> > boundary_vertices_idx;

 // get the largest boundary vertices as a ordered list (so you also have border edges)
 igl::boundary_loop(F, boundary_vertices_idx);

 // loop over each boundary loop 
 // taking only the largest boundary loop doesnt make sense here as 
 // stoke's law applies (unlike boundary conditions for the case of Tutte).
 // eigenvalue finding doesnt coverge if only largest boundary is taken
 for (int loopnum = 0; loopnum<boundary_vertices_idx.size(); loopnum++) {
  for (int item = 0; item<boundary_vertices_idx[loopnum].size(); item++) {
    int edg_i = boundary_vertices_idx[loopnum][item];
    int edg_j = boundary_vertices_idx[loopnum][(item + 1) % boundary_vertices_idx[loopnum].size()];

    triplets.push_back(Eigen::Triplet<double>(edg_i, edg_j + V_size, 0.5));
    triplets.push_back(Eigen::Triplet<double>(edg_i + V_size, edg_j, -0.5));
    triplets.push_back(Eigen::Triplet<double>(edg_j, edg_i + V_size, -0.5));
    triplets.push_back(Eigen::Triplet<double>(edg_j + V_size, edg_i, 0.5));
  }
 }

 A.setFromTriplets(triplets.begin(), triplets.end());


}

