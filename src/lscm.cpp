#include "lscm.h"
#include "vector_area_matrix.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/repdiag.h"
#include "igl/eigs.h"
#include <Eigen/SVD>
#include "igl/vector_area_matrix.h"
#include <iostream>
//#include <Eigen/Core>
//#include <Eigen/SparseCore>
//#include <SymGEigsSolver.h>
//#include <MatOp/SparseCholesky.h>
//#include <MatOp/SparseSymMatProd.h>
//#include <MatOp/SparseRegularInverse.h>


using namespace Eigen;

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	int n = V.rows();

	SparseMatrix<double> A;
	vector_area_matrix(F, A);
	//igl::vector_area_matrix(F, A);

	SparseMatrix<double> Q(2 * n, 2 * n), L(n, n);

	igl::cotmatrix(V, F, L);

	igl::repdiag(L, 2, Q);
	// Acc to assignment writeup
	Q -= 2*A;
	// Acc to Michael's comment
	//Q += 2*A;

	SparseMatrix<double> M(n, n), B(2 * n, 2 * n);

	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	igl::repdiag(M, 2, B);

	MatrixXd eigVec;
	VectorXd eigVal;

	//std::cout << Q << B << std::endl;

	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eigVec, eigVal);

	//Spectra::SparseSymMatProd<double> opQ(Q);
	//Spectra::SparseCholesky<double> opB(B);

	//Spectra::SymGEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEIGS_CHOLESKY> solver(&opQ, &opB, 3, 6);

	//solver.init();
	//int nconv = solver.compute(1000, 1e-10, Spectra::SMALLEST_MAGN);

	//assert(solver.info() == Spectra::SUCCESSFUL || nconv >= 3);
	//	eigVec = solver.eigenvectors();
	//	eigVal = solver.eigenvalues();

	U.resize(n, 2);
	U.col(0) = eigVec.col(2).head(n);
	U.col(1) = eigVec.col(2).tail(n);

	JacobiSVD <MatrixXd> svd(U, ComputeFullV);

	U = U * svd.matrixV();
}
