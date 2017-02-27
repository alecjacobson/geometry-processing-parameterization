#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/eigs.h>
#include <igl/cotmatrix.h>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  /*   
    Perform the following minimization:

    min ½ ∫ ‖ ∇u - ∇v⟂ ‖² dA
    
    which can be expanded into 

    min ½ ∫ ‖∇u‖² + ‖∇v‖² - ∇u·∇v⟂ dA
    
    min ½ DirichletNRG(u, v) - ∫∇u·∇v⟂ dA   

    min ½ DirichletNRG(u, v) - UᵀAU

    where the A comes from the vector_area_matrix calculation. The Dirichlet energy is formed through the cotangent laplacian.
    The minimization can be re-written as:
    
    min Uᵀ( / L 0 \ - A )U
            \ 0 L /
    
    this minimization will find the trivial solutions with zero cost as when: u→0, v→C OR u→C, v→0,  C≠0
    so we need to add additional constraints:

    1) minimize the energy
          don't mess with the given energy as defined above

    2) the solution has a non-zero norm
          ensure that the U is not set to a single point:
	  ∫‖u‖²dA = 1
	  where in our discrete case involves the mass matrix:
	  Uᵀ / M 0 \ U = 1
	     \ 0 M /
	  which is a quadratic constraint, which can be easily dealt with in a generalized eigenvalue problem

    3) is orthagonal to trivial solutions
          is orthagonal to trivial solutions, with the two trivial solutions energies (and e-vals) λ₁ = λ₂ = 0. Take the next non-zero eigenvalue as the solution
   */

  int32_t n = V.rows();
  U = Eigen::MatrixXd(n, 2);
  
  // construct the vector_area matrix:
  Eigen::SparseMatrix<double> A;
  vector_area_matrix(F, A);
  
  // construct the cotangent laplacian matrix:
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> Q;
  igl::repdiag(L, 2, Q);
  Q = Q - A;

  // construct the constraint matrix:
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  Eigen::SparseMatrix<double> B;
  igl::repdiag(M, 2, B);

  // solve the generalized eigenvalue problem, but only need the top 3 eigenvalues (first 2 are 0, 3rd is of interest)
  // according to the igl::eigs documentation only the EIGS_TYPE_SM is well supported
  Eigen::MatrixXd EVec;
  Eigen::VectorXd EVal;
  igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, EVec, EVal);
  
  Eigen::VectorXd UV = EVec.col(2);
  U.col(0) = UV.topRows(n);
  U.col(1) = UV.bottomRows(n);

  // get the canonical rotation:
  // models with high reflective symmetry will have prominent direction associated. Line them up using SVD?
  Eigen::MatrixXd UTU = U.transpose()*U;
  Eigen::JacobiSVD<Eigen::Matrix2d> SVD;
  SVD.compute(UTU, Eigen::ComputeFullV);
  Eigen::MatrixXd rotation = SVD.matrixV();

  U = U*rotation;
}
