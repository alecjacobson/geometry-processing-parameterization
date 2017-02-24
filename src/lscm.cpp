#include "lscm.h"
#include "vector_area_matrix.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/repdiag.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/eigs.h"

using namespace Eigen;

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	int n = V.rows();

	SparseMatrix<double> A;
	vector_area_matrix(F, A);

	SparseMatrix<double> Q(2 * n, 2 * n), L(n, n);

	igl::cotmatrix(V, F, L);

	igl::repdiag(L, 2, Q);
	// Acc to assignment writeup
	//Q -= A;
	// Acc to Michael's comment
	Q += A;

	SparseMatrix<double> M(n, n), B(2 * n, 2 * n);

	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);

	igl::repdiag(M, 2, B);

	MatrixXd eigVec;
	VectorXd eigVal;

	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eigVec, eigVal);

	VectorXd v = eigVec.col(2);

	U.resize(n, 2);
	U.col(0) = v.head(n);
	U.col(1) = v.tail(n);
}
