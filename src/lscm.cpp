#include "lscm.h"
#include "vector_area_matrix.h"

#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/repdiag.h"
#include "igl/eigs.h"
#include <Eigen/SVD>

using namespace std;

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int num_nodes = V.rows();
  U.resize(num_nodes, 2);

  Eigen::SparseMatrix<double> SignedArea;
  vector_area_matrix(F, SignedArea);

  // Obtain the cotangent Laplacian which covers the Dirichlet Energy terms
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> L_diag;
  igl::repdiag(L, 2, L_diag);
  // We factorize the 1/2 out to follow the Generalized Eigenvalue problem convention
  Eigen::SparseMatrix<double> Q = L_diag - 2 * SignedArea;

  // Mass Matrix used to enforce unit norm in our solution (this is our constraint)
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  Eigen::SparseMatrix<double> B;
  igl::repdiag(M, 2, B);
  
  // Borrowed from the libigl tutorial "Eigen Decomposition" section
  Eigen::MatrixXd EigenVectors;
  Eigen::VectorXd EigenValues;
  igl::eigs(Q, B, V.cols(), igl::EIGS_TYPE_SM, EigenVectors, EigenValues);

  Eigen::VectorXd FiedlerVector = EigenVectors.col(2);

  // Set U to the non-trivial eigen vector    
  U.col(0) = FiedlerVector.head(num_nodes);
  U.col(1) = FiedlerVector.tail(num_nodes);

  // Perform SVD to re-orient (PCA rotation is U)
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(U, Eigen::ComputeThinU | Eigen::ComputeThinV);
  U = svd.matrixU();
}
