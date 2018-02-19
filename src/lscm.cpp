#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include <igl/repdiag.h>
using namespace std;

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = V.leftCols(2);

  // My code
  int numV = V.rows();

  // construct B (mass) matrix
  Eigen::SparseMatrix<double> M;
	igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  Eigen::SparseMatrix<double> B(V.rows()*2, V.rows()*2);
  igl::repdiag(M,2,B);

  // construct A (vector area) matrix
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F,A);

  // construct LL (concatenated laplace) matrix
  Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V,F,L);
  Eigen::SparseMatrix<double> LL(V.rows()*2, V.rows()*2);
  igl::repdiag(L,2,LL);

  // construct Q (LL - A) matrix
  Eigen::SparseMatrix<double> Q;
  Q = LL - A;

  // eigs
  Eigen::MatrixXd eVec;
	Eigen::VectorXd eVal;
	// cout << "solving eigs ... (may takes a while)" << endl;
  igl::eigs(Q, B, 5, igl::EIGS_TYPE_SM, eVec, eVal);
 	// cout << "finish eigs" << endl;

 	// assign to U
 	U.resize(numV, 2);
 	for (int ii = 0; ii<numV; ii++){
 		U(ii,0) = eVec(ii,1);
 		U(ii,1) = eVec(ii+numV,1);
 	}
}
