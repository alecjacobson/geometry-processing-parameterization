#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

void vector_area_matrix(
        const Eigen::MatrixXi &F,
        Eigen::SparseMatrix<double> &A) {

    int n = F.maxCoeff() + 1;
    A.resize(n * 2, n * 2);

    // In this case, we want a list of boundaries
    // where L[i] = ordered list of boundary vertices in loop i
    // Given this is an overloaded function, let's pass a std::vector as second argument.
    std::vector<std::vector<int>> bnds;
    igl::boundary_loop(F, bnds);

    typedef Eigen::Triplet<double> eT;
    std::vector<eT> list;

    for (auto bnd : bnds) {
        for (int j = 0; j < bnd.size(); j++) {
            int u_i = bnd[j];
            int u_j = bnd[(j + 1) % bnd.size()];
            list.emplace_back(eT(u_i, u_j + n, 0.5));
            list.emplace_back(eT(u_i + n, u_j, -0.5));
        }
    }
    Eigen::SparseMatrix<double> A_(n * 2, n * 2);
    A_.setFromTriplets(list.begin(), list.end());

    Eigen::SparseMatrix<double> A_T(n * 2, n * 2);
    A_T = A_.transpose();
    A = (A_ + A_T);

}

