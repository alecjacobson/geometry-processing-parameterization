#include "lscm.h"
#include "vector_area_matrix.h"
#include "igl/eigs.h"
#include "igl/repdiag.h"
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"

#include "igl/lscm.h"

void lscm(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &U) {

    int n = V.rows();
    U.resize(n, 2);

    Eigen::SparseMatrix<double> A;
    vector_area_matrix(F, A);

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    Eigen::SparseMatrix<double> Q;
    igl::repdiag(L, 2, Q);
    Q = Q - A;


    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    Eigen::SparseMatrix<double> B;
    igl::repdiag(M, 2, B);

    Eigen::MatrixXd eigen_vectors;
    Eigen::VectorXd eigen_values;

    // The two firsts are trivial solutions.
    // so, let's use k=3 and get the the last one.
    int k = 3;
    igl::eigs(Q, B, k, igl::EIGS_TYPE_SM, eigen_vectors, eigen_values);

    U.col(0) = eigen_vectors.col(k-1).head(n);
    U.col(1) = eigen_vectors.col(k-1).tail(n);

    // rotation
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose()*U, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = U * svd.matrixU();


}
