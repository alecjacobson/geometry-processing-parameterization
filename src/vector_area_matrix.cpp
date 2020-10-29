#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
    int V_size = F.maxCoeff() + 1;
    A.resize(V_size * 2, V_size * 2);

    std::vector<std::vector<int>> b;
    igl::boundary_loop(F, b);

    std::vector<Eigen::Triplet<double>> tripletList;

    for (int i = 0; i < b.size(); i++) {
        std::vector<int> loop = b[i];

        for (int j = 0; j < loop.size(); j++) {

            int ui = loop[j];
            int uj = loop[(j + 1) % loop.size()];

            tripletList.push_back(Eigen::Triplet<double>(ui + V_size, uj, -0.25));
            tripletList.push_back(Eigen::Triplet<double>(uj, ui + V_size, -0.25));
            tripletList.push_back(Eigen::Triplet<double>(ui, uj + V_size, 0.25));
            tripletList.push_back(Eigen::Triplet<double>(uj + V_size, ui, 0.25));
        }
    }

    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

