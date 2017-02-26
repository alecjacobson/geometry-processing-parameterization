#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

using namespace Eigen;
using namespace std;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
	// Replace with your code

	int n = F.maxCoeff() + 1;
	vector<Triplet<double>> coef;
	vector<vector<int>> bndr;
	igl::boundary_loop(F, bndr);

	for (auto& b : bndr) 
	{
		int bsize = b.size();
		for (int v = 0; v < bsize; ++v)
		{
			int v1 = b[v], v2 = b[(v + 1) % bsize];
	
			coef.push_back(Triplet<double>(v1, v2 + n, 0.25));
			coef.push_back(Triplet<double>(v2, v1 + n, -0.25));
			coef.push_back(Triplet<double>(v2 + n, v1, 0.25));
			coef.push_back(Triplet<double>(v1 + n, v2, -0.25));			
		}		
	}
	
	A.resize(2 * n, 2 * n);
	A.setZero();
	A.setFromTriplets(coef.begin(), coef.end());
}

