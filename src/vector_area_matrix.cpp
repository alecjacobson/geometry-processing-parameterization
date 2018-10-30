#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int V_size = F.maxCoeff()+1;
  A.resize(V_size*2,V_size*2);

  // Determine the longest boundary loop
  Eigen::VectorXi bdry;
  igl::boundary_loop(F, bdry);
  
  // Each edge (i, j) in the loop contributes 1/2 (x_i y_j - x_j y_i) to the vector area
  typedef Eigen::Triplet<double> T;
  std::vector<T> triples;
  triples.reserve(bdry.size() * 4);
  for (int i = 0; i < bdry.size(); i++) {
	if (i == bdry.size() - 1) {
		triples.push_back(T(bdry[i], bdry[0] + V_size, 0.25));
		triples.push_back(T(bdry[0] + V_size, bdry[i], 0.25));
   		triples.push_back(T(bdry[0], bdry[i] + V_size, -0.25));
   		triples.push_back(T(bdry[i] + V_size, bdry[0], -0.25));
	} 
	else {
		triples.push_back(T(bdry[i], bdry[i+1] + V_size, 0.25));
		triples.push_back(T(bdry[i+1] + V_size, bdry[i], 0.25));
   		triples.push_back(T(bdry[i+1], bdry[i] + V_size, -0.25));
   		triples.push_back(T(bdry[i] + V_size, bdry[i+1], -0.25));
	}
  }
  A.setFromTriplets(triples.begin(), triples.end());
}

