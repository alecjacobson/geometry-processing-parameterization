#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
	// Compute cotangent laplacian
	typedef Eigen::SparseMatrix<double> SpMat;
	SpMat L;
	igl::cotmatrix(V, F, L);

	// Compute Q in the U^TQU for the minimization
	//		L 0
	// Q = (   ) - A
	//      0 L 
	SpMat repL, A;
	igl::repdiag(L, 2, repL);
	vector_area_matrix(F, A);
	SpMat Q = repL - 2 * A;

	// Compute mass matrix and from B in constraint U^TBU = 1
	//		M 0
	// B = (   ) where M is mass matrix
	//      0 M 
	SpMat M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	SpMat B;
	igl::repdiag(M, 2, B);

	// Implement eigen-decomposition to find the eigenvector with the 
	// 3rd smallest eigenvalue which satisfies all criteria
	Eigen::MatrixXd eigenvector;
	Eigen::VectorXd eigenvalue;
	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eigenvector, eigenvalue);
	int n = F.maxCoeff() + 1;
	U.resize(n, 2);
	U.col(0) = eigenvector.col(2).head(n);
	U.col(1) = eigenvector.col(2).tail(n);

	// Find canonical rotation using singular value decomposition
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose()*U, Eigen::ComputeThinU | Eigen::ComputeThinV);
	U = U * svd.matrixU();
}
