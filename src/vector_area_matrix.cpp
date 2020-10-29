#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
	const Eigen::MatrixXi & F,
	Eigen::SparseMatrix<double>& A)
{
	int V_size = F.maxCoeff() + 1;

	std::vector<std::vector<int>> boundary;
	igl::boundary_loop(F, boundary);
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	for (int b = 0; b < boundary.size(); b++) {
		std::vector<int> bound = boundary[b];
		for (int i = 0; i < bound.size(); i++) {
			int j = (i + 1) % bound.size();
			tripletList.push_back(T(bound[i], bound[j] + V_size, 0.25));
			tripletList.push_back(T(bound[i] + V_size, bound[j], -0.25));
			tripletList.push_back(T(bound[j], bound[i] + V_size, -0.25));
			tripletList.push_back(T(bound[j] + V_size, bound[i], 0.25));
		}
	}
	A.resize(V_size * 2, V_size * 2);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
}

