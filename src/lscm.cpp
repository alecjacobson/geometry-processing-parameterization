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
	//U = V.leftCols(2);
	int n = V.rows();

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	Eigen::SparseMatrix<double> VA;
	//vector_area_matrix(F, VA);
	igl::vector_area_matrix(F, VA);

	Eigen::SparseMatrix<double> L_rep;
	igl::repdiag(L, 2, L_rep);

	Eigen::SparseMatrix<double> Q = L_rep - 2*VA;

	Eigen::SparseMatrix<double> B;
	igl::repdiag(M, 2, B);

	Eigen::MatrixXd E_vec;
	Eigen::VectorXd E_val;
	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, E_vec, E_val);

	std::cout << E_val << std::endl;
	Eigen::VectorXd vectorized_U = E_vec.col(2);
	U.col(0) = vectorized_U.topRows(n);
	U.col(1) = vectorized_U.bottomRows(n);
}
