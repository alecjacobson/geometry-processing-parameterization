#include "lscm.h"

#include <vector>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include "vector_area_matrix.h"
#include <igl/vector_area_matrix.h>
#include <igl/repdiag.h>

#include <igl/lscm.h>
#include <igl/boundary_loop.h>
void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	// Replace with your code
	int n = V.rows();
	
	Eigen::SparseMatrix<double> A, L(n,n), Q(2*n, 2*n), M(n,n), B(2*n,2*n);

	vector_area_matrix(F, A);	//igl::vector_area_matrix(F, A); 
	igl::cotmatrix(V, F, L);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	
	igl::repdiag(L, 2, Q);
	
	//Convert A to a symetric matrix and apply it to Q. Our implementation
	//is OK, no need to convert.
	//Eigen::SparseMatrix<double, Eigen::ColMajor> A2(A.transpose());
	Q -= A;	
	igl::repdiag(M, 2, B);
	
	Eigen::MatrixXd eigenVector;
	Eigen::VectorXd eigenValue;		
	
	igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, eigenVector, eigenValue);

	//The first 2 eigen vectors aren't the best solutions, and still give 
	//overlapping regions in UV space. The 3rd eigen vector gives the required
	//results

	U.resize(n, 2);
	U.col(0) = eigenVector.col(2).head(n);
	U.col(1) = eigenVector.col(2).tail(n);
	
	//Do this cannonical rotation. Not sure if this is making a difference,
	//in my parameterization.
	Eigen::JacobiSVD <Eigen::MatrixXd> svd(U, Eigen::ComputeFullV);	
	U = U * svd.matrixV();

	////THIS WORKS, IF ALL ELSE FAILS
	//Eigen::VectorXi bnd, b(2, 1);
	//igl::boundary_loop(F, bnd);
	//b(0) = bnd(0);

	//b(1) = bnd(round(bnd.size() / 2));
	//Eigen::MatrixXd bc(2, 2);
	//bc << 0, 0, 1, 0;

	//igl::lscm(V, F, b, bc, U);
}




