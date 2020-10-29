#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

typedef Eigen::Triplet<double> T;

void vector_area_matrix(const Eigen::MatrixXi & F,
		Eigen::SparseMatrix<double>& A) {
	// Replace with your code
	int V_size = F.maxCoeff() + 1;
	A.resize(V_size * 2, V_size * 2);

	Eigen::VectorXd bounds;
	igl::boundary_loop(F, bounds);
	int bnum = bounds.size();

	std::vector<T> coef;
	for (int i = 0; i < bnum; i++) {
		int b1 = bounds[i];
		for (int j = 0; j < bnum; j++) {
			int b2 = bounds[j];
			T tmp(b1, b2 + V_size, 0.5);
			coef.push_back(tmp);

			T tmp2(b1 + V_size, b2, -0.5);
			coef.push_back(tmp2);
		}
	}

	Eigen::SparseMatrix<double> A_(V_size* 2, V_size*2);
	A_.setFromTriplets(coef.begin(), coef.end());
	A = A_ + A_.transpose();
}

