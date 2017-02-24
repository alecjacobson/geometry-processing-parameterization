#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
void vector_area_matrix(
	const Eigen::MatrixXi & F,
	Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
	const auto v = F.maxCoeff() + 1;
	A.resize(2*v, 2*v);

	std::vector<Eigen::Triplet<double>> entries;

	std::vector<std::vector<int>> boundaries;
	igl::boundary_loop(F, boundaries);
	for (const auto& loop : boundaries) {
		const int loop_length = loop.size();
		for (int vert = 0; vert < loop_length; vert++) {
			//Have an edge from loop[i], loop[i+1 % loop.size()]
			int i = loop[vert];
			int j = loop[(vert + 1) % loop_length];

			entries.emplace_back(i, j + v, 0.25);
			entries.emplace_back(j + v, i, 0.25);
			entries.emplace_back(i + v, j, -0.25);
			entries.emplace_back(j, i + v, -0.25);
		}
	}

	A.setFromTriplets(entries.begin(), entries.end());
}

