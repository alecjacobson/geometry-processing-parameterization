#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  //A.resize(F.maxCoeff()+1,F.maxCoeff()+1);

	int n = F.maxCoeff() + 1;
	A.resize(2 * n, 2 * n);	
	
	std::vector<std::vector<int>> boundaries;
	igl::boundary_loop(F, boundaries);
	std::vector<Eigen::Triplet<double>> triplets;

	for (int i = 0; i < boundaries.size(); i++) 
	{
		std::vector<int> bound = boundaries[i];
		for (int j = 0; j < bound.size();j++)
		{
			int x = bound[j];
			int y = bound[(j + 1) % bound.size()];

			//Create A |ui uj| -> uix*ujy - ujx*uiy
			triplets.emplace_back(y + n, x, 1);
			triplets.emplace_back(x + n, y, -1);

			//Create A^T
			triplets.emplace_back(x, y + n, 1);
			triplets.emplace_back(y, x + n, -1);
		}
	}	
	
	A.setFromTriplets(triplets.begin(), triplets.end());
	A = 0.5*A; // 1/2(A + A^T)
}

