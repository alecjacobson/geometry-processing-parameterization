#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

using namespace Eigen;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
	std::vector<std::vector<int>> boundaries;
	igl::boundary_loop(F, boundaries);

	std::vector<Triplet<double>> A_val;

	auto reservation = [&boundaries]()->size_t {
		size_t res = 0;
		for each(auto b in boundaries)
			res += b.size();

		return 4 * res;
	};

	A_val.reserve(reservation());

	int n = F.maxCoeff() + 1;

	for each (auto b in boundaries)
	{
		int v1 = b[0], v2;
		for (int i = 1; i < b.size(); ++i)
		{
			v2 = b[i];

			A_val.push_back({ v1, n + v2, 0.5 });
			A_val.push_back({ v2, n + v1, -0.5 });
			A_val.push_back({ n + v2, v1, 0.5 });
			A_val.push_back({ n + v1, v2, -0.5 });

			v1 = v2;
		}
		v2 = b[b.size() - 1];
	}

	A.resize(2 * n, 2 * n);
	A.setZero();
	A.setFromTriplets(A_val.begin(), A_val.end());
}

