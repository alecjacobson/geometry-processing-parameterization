#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
#include <iostream>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  int num_vect = F.maxCoeff() + 1;
  std::vector<std::vector<int>> L;
  igl::boundary_loop(F, L);
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  for (int i = 0; i < L.size(); i++){
	for (int j = 0; j < L[i].size(); j++){ 
	 	int a = L[i][j];
	 	int k = (j + 1) % L[i].size();
	 	int b = L[i][k];

	 	tripletList.push_back(T(a, b + num_vect, 1));
        tripletList.push_back(T(a + num_vect, b, -1));
	 	}
	 }

  A.resize(2 * num_vect, 2 * num_vect);
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  Eigen::SparseMatrix<double> AT = A.transpose();
  A = .5 * (A + AT);
  std::cout << "here" << std::endl;
}