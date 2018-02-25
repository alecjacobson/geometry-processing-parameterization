#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
	// Compute all boundary edges
	std::vector<std::vector<int>> L;
	igl::boundary_loop(F, L);

	// Construct A so that ½∑_{{i,j}∈∂S}|ui uj| can be written as U^TAU with U = [u;v]
	int n = F.maxCoeff() + 1;
	A.resize(n * 2, n * 2);
	typedef Eigen::Triplet<double> T;
	typedef Eigen::SparseMatrix<double> SpMat;
	std::vector<T> triplets;
	for (int i = 0; i < L.size(); i++) {
		for (int j = 0; j < L[i].size(); j++) {
			int u_i = L[i][j];
			int u_j;
			if (j == (L[i].size() - 1))
				u_j = L[i][0];
			else
				u_j = L[i][j + 1];
			triplets.push_back(T(u_i, u_j + n, 0.5));
			triplets.push_back(T(u_i + n, u_j, -0.5));
		}
	}
	// Make sure A is symmetric
	SpMat A_hat(n * 2, n * 2), A_hatT(n * 2, n * 2);
	A_hat.setFromTriplets(triplets.begin(), triplets.end());
	A_hatT = A_hat.transpose();
	A = 0.5*(A_hat + A_hatT);
}

