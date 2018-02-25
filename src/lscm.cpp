#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>
#include <Eigen/SVD>
#include <Eigen/LU>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // // Replace with your code
  // U = V.leftCols(2);


	// We need to make the matrices that correspond to matrices Q and B in the readme.
	// First, make Q by computing the cotangent matrix.

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	// Now use igl's repdiag function to make it twice as large along the diagonal
	int N_verts = F.maxCoeff()+1;
	Eigen::SparseMatrix<double> L_big (N_verts*2, N_verts*2);
	igl::repdiag(L, 2, L_big);

	// Create vector area matrix
	Eigen::SparseMatrix<double> A;
  	vector_area_matrix(F, A);

	// Assemble Q
	Eigen::SparseMatrix<double> Q (N_verts*2, N_verts*2);
	Q = L_big - A;

	// Now work on B. First, get the mass matrix. Use the default type.
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	// Make its bigger version, B
	Eigen::SparseMatrix<double> B (N_verts*2, N_verts*2);
	igl::repdiag(M, 2, B);

	// Now we have everything we need to find the eigensolutions of the system as per the readme.
	// Based on the readme, we want to ignore the 2 smallest (~0) eigenvalues, and get the third one.
	Eigen::MatrixXd sU;
	Eigen::VectorXd sS;
	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, sU, sS);

	// Extract the third eigenvector
	// Eigen::MatrixXd eig_vec_3 = sU.col(2);

	// Put these values in terms of u and v coordinates into the matrix U
	U.resize(N_verts, 2);
	U.col(0) = sU.col(2).block(0, 0, N_verts, 1);
	U.col(1) = sU.col(2).block(N_verts, 0, N_verts, 1);

	// To do the SVD, compute U^T*U, which is the covariance matrix
	Eigen::MatrixXd UTU = U.transpose()*U;

	// Compute the SVD of this matrix
	Eigen::JacobiSVD<Eigen::Matrix2d> svd(UTU, Eigen::ComputeFullU | Eigen::ComputeFullV);

	// Now the V matrix should provide the rotation we need, to compute the rotated matrix
	U = U*svd.matrixV();

	return;

}
  


