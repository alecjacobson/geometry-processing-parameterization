#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
#include <vector>
void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
	int n = F.maxCoeff() + 1;
	
	// Replace with your code
	std::vector<std::vector<int>> boundaries;
	igl::boundary_loop(F, boundaries);
	int bn = boundaries.size();

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(4 * bn);
	// Build 1/2(A + A^T) to make it symmetric. 
	// Loop around for last element
	std::vector<int> boundary;
	for (int k = 0; k < bn; ++k) 
	{
		boundary = boundaries[k];
		int length = boundary.size();
		for (int l = 0; l < length; ++l) 
		{
			int i = boundary[l];
			int j = boundary[(l + 1) % length];

			//A
			triplets.push_back(Eigen::Triplet<double>(j + n, i, 1));
			triplets.push_back(Eigen::Triplet<double>(i + n, j, -1));

			//A^T
			triplets.push_back(Eigen::Triplet<double>(i, j + n, 1));
			triplets.push_back(Eigen::Triplet<double>(j, i + n, -1));
		}
		
	}

	A.setFromTriplets(triplets.begin(),triplets.end());
	A = 0.5*A;


}

