#include "lscm.h"
#include "vector_area_matrix.h"
#include "igl/massmatrix.h"
#include "igl/repdiag.h"
#include "igl/eigs.h"
#include "igl/cotmatrix.h"

using namespace Eigen;
using namespace igl;

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    // Replace with your code
	int nvert = V.rows();

	// Computing Q matrix
	
	SparseMatrix<double> L;
	cotmatrix(V, F, L);

	SparseMatrix<double> KK(2 * nvert, 2 * nvert);
	repdiag(L, 2, KK);
	
	SparseMatrix<double> A;
	vector_area_matrix(F, A);
	SparseMatrix<double> Q = KK - 2*A;
	
	// Computing constraint matrix

	SparseMatrix<double> M(nvert, nvert);
	massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	SparseMatrix<double> B(2 * nvert, 2 * nvert);
	repdiag(M, 2, B);

	// Solving Generalized Eigenvalue problem

	int k = 4;
	MatrixXd sU;
	VectorXd sS;
	eigs(Q, B, k, igl::EIGS_TYPE_SM, sU, sS);

	// Find non zero eigen value

	double nonZeroValue;
	for (int i = sS.size() - 1; i >= 0; --i)
	{
		if (sS[i] > 0) 
		{
			nonZeroValue = sS[i];
		}
	}
	
}
