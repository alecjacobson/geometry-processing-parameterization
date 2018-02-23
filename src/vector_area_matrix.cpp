#include "vector_area_matrix.h"

#include "igl/boundary_loop.h"

#include <vector>

typedef Eigen::Triplet<double> T;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int n = F.maxCoeff() + 1;
  A.resize(2 * n, 2 * n);

  std::vector<std::vector<int>> boundary_loops;
  igl::boundary_loop(F, boundary_loops);

	std::vector<T> tripletList;
	for(int i = 0; i < boundary_loops.size(); i++)
	{
		for(int j = 0; j < boundary_loops[i].size(); j++) {
			int ui = boundary_loops[i][j];
			int vi = boundary_loops[i][j] + n;
			int uj = boundary_loops[i][(j + 1) % boundary_loops[i].size()];
			int vj = boundary_loops[i][(j + 1) % boundary_loops[i].size()] + n;

			tripletList.push_back(T(ui, vj, 0.5));	
			tripletList.push_back(T(uj, vi, -0.5));	
		}
	}

	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A += Eigen::SparseMatrix<double>(A.transpose());
}

