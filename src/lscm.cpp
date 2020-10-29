#include "lscm.h"
#include "vector_area_matrix.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/repdiag.h"
#include "igl/eigs.h"

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    Eigen::SparseMatrix<double> A;
    vector_area_matrix(F, A);

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    Eigen::SparseMatrix<double> L_diag;
    igl::repdiag(L, 2, L_diag);

    Eigen::SparseMatrix<double> Q = L_diag - 2. * A;

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    Eigen::SparseMatrix<double> B;
    igl::repdiag(M, 2, B);

    Eigen::MatrixXd v;
    Eigen::VectorXd lambda;
    igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, v, lambda);

    const int size = V.rows();
    U.resize(size, 2);
    U.col(0) = v.col(2).head(size);
    U.col(1) = v.col(2).tail(size);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose() * U, Eigen::ComputeThinU | Eigen::ComputeThinV);

    U = U * svd.matrixU();
}
