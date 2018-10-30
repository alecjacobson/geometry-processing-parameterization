#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int V_size = F.maxCoeff()+1;
  Eigen::SparseMatrix<double> A_u;
  A_u.resize(V_size*2,V_size*2);

  // Determine the longest boundary loop
  std::vector<std::vector<int>> bdry_loops;
  igl::boundary_loop(F, bdry_loops);
  
  // Each edge (i, j) in the loop contributes 1/2 (x_i y_j - x_j y_i) to the vector area
  typedef Eigen::Triplet<double> T;
  std::vector<T> triples;
  triples.reserve(bdry_loops[0].size() * 4);
  for (int l = 0; l < bdry_loops.size() ; l++) {
  std::vector<int> bdry = bdry_loops[l];
  for (int i = 0; i < bdry.size(); i++) {
	if (i == bdry.size() - 1) {
		triples.push_back(T(bdry[i], bdry[0] + V_size, 0.5));
   		triples.push_back(T(bdry[i] + V_size, bdry[0], -0.5));
	} 
	else {
		triples.push_back(T(bdry[i], bdry[i+1] + V_size, 0.5));
   		triples.push_back(T(bdry[i] + V_size, bdry[i+1], -0.5));
	}
  }
}
  A_u.setFromTriplets(triples.begin(), triples.end());
  Eigen::SparseMatrix<double> B;
  B = A_u.transpose();
  A = (A_u + B) / 2;
}

