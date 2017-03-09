#include "lscm.h"
#include <vector_area_matrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/vector_area_matrix.h>

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

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	Eigen::SparseMatrix<double> VA(2 * n,2 * n);
	vector_area_matrix(F, VA);
	//igl::vector_area_matrix(F, VA);

	Eigen::SparseMatrix<double> L_rep(2 * n, 2 * n);
	igl::repdiag(L, 2, L_rep);

	Eigen::SparseMatrix<double> Q = L_rep - VA;

	Eigen::SparseMatrix<double> B(2 * n, 2 * n);
	igl::repdiag(M, 2, B);

	// Solve generalized eigenvalue problem:
	// Note: this requires a change to eigs.cpp that won't be submitted alongside my code due to .gitignore
	// See here https://github.com/alecjacobson/geometry-processing-parameterization/issues/3
	Eigen::MatrixXd E_vec;
	Eigen::VectorXd E_val;
	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, E_vec, E_val);

	// Select first non-zero eigenvalue
	Eigen::VectorXd vectorized_U = E_vec.col(2);
	U.col(0) = vectorized_U.topRows(n);
	U.col(1) = vectorized_U.bottomRows(n);

	// Rotate parameterized mesh by principle direction
	Eigen::MatrixXd C = U.transpose()*U;
	Eigen::JacobiSVD<Eigen::Matrix2d> svd;
	svd.compute(C, Eigen::ComputeFullV);
	Eigen::MatrixXd W = svd.matrixV();
	U = U*W;

}
