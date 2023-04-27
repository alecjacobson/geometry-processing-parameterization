#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"
#include <iostream>

using namespace std;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  int V_size = F.maxCoeff()+1;
  A.resize(V_size*2,V_size*2);
  A.setZero();
  Eigen::SparseMatrix<double> A_bar(V_size * 2, V_size * 2);
  std::vector<Eigen::Triplet<double>> triplets;
  // Detect the boundary vertices
  std::vector<std::vector<int>> lb;
  igl::boundary_loop(F, lb);

  for (std::vector<std::vector<int>>::iterator it = lb.begin(); it != lb.end(); ++it) {
	  for (vector<int>::iterator bdry = (*it).begin(); bdry != (*it).end(); ++bdry) {
		  int i, j;
		  if ((bdry + 1) != (*it).end()) {
			  i = *bdry; j = *(bdry+1);
		  } else {
			i = *bdry; j = *((*it).begin());
		  }
		  // Fill in entries of A_bar as well as entries of A_bar ^ T
		  triplets.push_back({ i, j + V_size,  1./2. });
		  triplets.push_back({ i + V_size, j, -1./2. });
		  triplets.push_back({ j + V_size, i,  1./2. });
		  triplets.push_back({ j, i + V_size, -1./2. });
	  }
  }
  A.setFromTriplets(triplets.begin(), triplets.end());
}
