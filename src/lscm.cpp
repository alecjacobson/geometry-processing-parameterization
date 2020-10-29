#include "lscm.h"
#include <igl/cotmatrix.h>
#include "vector_area_matrix.h"
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    Eigen::SparseMatrix<double> L, A, Q, M, B;
    Eigen::MatrixXd evecs;
    Eigen::VectorXd evals;
    
    igl::cotmatrix(V, F, L);
    vector_area_matrix(F, A);
    igl::repdiag(L, 2, Q);
    Q -= A * 2;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::repdiag(M, 2, B);
    igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, evecs, evals);
    U.resize(V.rows(), 2);
    U.col(0) = evecs.col(2).head(V.rows());
    U.col(1) = evecs.col(2).tail(V.rows());
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose() * U, Eigen::ComputeThinU | Eigen::ComputeThinV);
    U *= svd.matrixU();
    U.col(0) = -U.col(0);
}
