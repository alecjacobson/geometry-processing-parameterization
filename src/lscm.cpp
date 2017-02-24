#include "lscm.h"
#include <vector>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include "vector_area_matrix.h"

void lscm(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::MatrixXd & U)
{
  // Replace with your code
	U = V.leftCols(2);
	int n = V.rows();

	Eigen::SparseMatrix<double> Q(2 * n, 2 * n);
	{
		Eigen::SparseMatrix<double> L;
		igl::cotmatrix(V, F, L);

		std::vector<Eigen::Triplet<double>> Q_entries;
		for (int k = 0; k < L.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
				Q_entries.emplace_back(it.row(), it.col(), it.value());
				Q_entries.emplace_back(it.row() + n, it.col() + n, it.value());
			}
		}
		Q.setFromTriplets(Q_entries.begin(), Q_entries.end());
	}
	Eigen::SparseMatrix<double> A;
	vector_area_matrix(F, A);
	Q += A;

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	Eigen::SparseMatrix<double> B(2 * n, 2 * n);
	
	std::vector<Eigen::Triplet<double>> B_entries;
	B_entries.reserve(2*M.nonZeros());
	B.reserve(Eigen::VectorXi::Constant(2 * n, 1));
	for (int i = 0; i < n;i++) {
		B.insert(i, i) = M.diagonal()(i);
		B.insert(i + n, i + n) = M.diagonal()(i);
	}
	B.makeCompressed();


	Eigen::VectorXd eigenvals;
	Eigen::MatrixXd eigenvecs;
	int k = 3;//2*n;


	igl::eigs(Q, B, k, igl::EIGS_TYPE_SM, eigenvecs, eigenvals);

	//Pick the smallest non-zero eigenvalue. They should be sorted in descending order, so start at the end.
	for (int i = k - 1; i >= 0; i--) {
		if (eigenvals(i) != 0) {
			U.col(0) = eigenvecs.block(0, i, n, 1);
			U.col(1) = eigenvecs.block(n, i, n, 1);
			break;
		}
	}

	//Finally, find the canonical rotation.
}
