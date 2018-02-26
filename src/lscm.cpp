#include "lscm.h"
#include <Eigen/Core>
#include "igl/cotmatrix.h"
#include "igl/repdiag.h"
#include "igl/massmatrix.h"
#include "igl/eigs.h"
#include "vector_area_matrix.h"

void lscm(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		Eigen::MatrixXd & U) {

	int vertexCount = F.maxCoeff() + 1;

	Eigen::SparseMatrix<double> areaMatrix;
	vector_area_matrix(F, areaMatrix);

	Eigen::SparseMatrix<double> lagrangian;
	igl::cotmatrix(V, F, lagrangian);

	// This is the symmetric quadratic matrix from the readme.
	Eigen::SparseMatrix<double> Q;
	igl::repdiag(lagrangian, 2, Q);
	Q = Q - areaMatrix;

	// Now we create the B matrix constraint from the readme.
	Eigen::SparseMatrix<double> massMatrix;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, massMatrix);
	Eigen::SparseMatrix<double> B;
	igl::repdiag(massMatrix, 2, B);

	// Now we solve using eigenvalue decomposition. We skip the first two
	// solution/eigenvector, thus we use k = 3.
	Eigen::MatrixXd eigenVectors;
	Eigen::VectorXd eigenValues;
	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eigenVectors, eigenValues);

	// Now we handle that problem with the rotational ambiguity by using
	// SVD. We take the third solution/eigenvector and split them by axes.
	Eigen::MatrixXd rawU(vertexCount, 2);
	rawU.col(0) = eigenVectors.col(2).head(vertexCount);
	rawU.col(1) = eigenVectors.col(2).tail(vertexCount);

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(rawU, Eigen::ComputeThinU | Eigen::ComputeThinV);
	U = svd.matrixU() * svd.singularValues().asDiagonal();
}
