#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>

void lscm(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & U)
{
    Eigen::SparseMatrix<double> L, A, Q, M, B;
    igl::cotmatrix(V, F, L);
    vector_area_matrix(F, A);

    igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
    Q = igl::repdiag(L, 2) - A;
    B = igl::repdiag(M, 2);

    Eigen::MatrixXd sU;
    Eigen::VectorXd sS;
    igl::eigs(Q, B, 3, igl::EIGS_TYPE_SM, sU, sS);

    U.resize(sU.rows() / 2, 2);
    Eigen::VectorXd u = sU.col(2) / sU.col(2).maxCoeff();
    U.col(0) = u.head(sU.rows() / 2);
    U.col(1) = u.tail(sU.rows() / 2);
}
