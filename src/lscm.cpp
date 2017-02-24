#include "lscm.h"
#include <vector>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include "vector_area_matrix.h"
#include <igl/vector_area_matrix.h>
#include <igl/repdiag.h>
void lscm(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::MatrixXd & U)
{
  // Replace with your code
	U = V.leftCols(2);
	int n = V.rows();

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	auto LD = igl::repdiag(L, 2);

	Eigen::SparseMatrix<double> A;
	vector_area_matrix(F, A);
	//igl::vector_area_matrix(F, A);

	Eigen::SparseMatrix<double> Q = LD + 2*A;

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	Eigen::SparseMatrix<double> B(2 * n, 2 * n);
	
	B = igl::repdiag(M, 2);

	Eigen::VectorXd eigenvals;
	Eigen::MatrixXd eigenvecs;
	int k = 3;


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
	//NOTE: Should center U first. 

	Eigen::RowVector2d mean = U.colwise().mean();
	U.rowwise() -= mean;

	Eigen::Matrix2d covariance = U.transpose()*U;
	Eigen::JacobiSVD<Eigen::Matrix2d> jacobi;
	jacobi.compute(covariance, Eigen::ComputeFullV);
	auto pca = jacobi.matrixV();
	//The columns of pca are the axes of the ideal 
	for (int i = 0; i < U.rows(); ++i) {
		U.row(i) = (pca*U.row(i).transpose()).transpose().eval();
	}
}
