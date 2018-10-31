#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>
#include <Eigen/SVD>

void lscm(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::MatrixXd & U)
{
	// construct Q and B
	Eigen::SparseMatrix<double> L, L2, A, Q, M, B;
	igl::cotmatrix(V, F, L);
	igl::repdiag(L, 2, L2);
	vector_area_matrix(F, A);
	Q.resize(A.rows(), A.rows()); // A is 2#V x 2#V
	Q = L2 - 2.0 * A;
	igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
	igl::repdiag(M, 2, B);

	// solve generalized Eigen value problem
	Eigen::MatrixXd eigenVec;
	Eigen::VectorXd eigenVal;
	igl::eigs(Q, B, 3, igl::EigsType::EIGS_TYPE_SM, eigenVec, eigenVal);
	
	U.resize(V.rows(), 2);
	U.col(0) = eigenVec.col(2).head(V.rows());
	U.col(1) = eigenVec.col(2).tail(V.rows());
	// rotate
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose() * U, Eigen::ComputeFullU | Eigen::ComputeFullV);
	U = U * svd.matrixU();
}
