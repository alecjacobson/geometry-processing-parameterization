#include "lscm.h"

#include "vector_area_matrix.h"

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>

void lscm(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,
		Eigen::MatrixXd & U) {
	// Replace with your code
	//  U = V.leftCols(2);

	// L
	// calculate L, here, w is 1
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	// A
	Eigen::SparseMatrix<double> A;
	vector_area_matrix(F, A);

	// Q
	Eigen::SparseMatrix<double> LL;
	LL = igl::repdiag(L, 2);
	Eigen::SparseMatrix<double> Q;
	Q = L - A;

	// B
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	Eigen::SparseMatrix<double> B;
	B = igl::repdiag(M, 2);

	// equation
	Eigen::MatrixXd Vecs;
	Eigen::VectorXd lambdas;
	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, Vecs, lambdas);

	int vnum = V.rows();
	U.resize(vnum, 2);
	U.col(0) = Vecs.col(2).head(vnum);
	U.col(1) = Vecs.col(2).tail(vnum);

	// final rotation
	Eigen::MatrixXd utu = U.transpose() * U;
	// svd
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(utu,
			Eigen::ComputeFullU | Eigen::ComputeFullV);
	U = U * svd.matrixU();
}
