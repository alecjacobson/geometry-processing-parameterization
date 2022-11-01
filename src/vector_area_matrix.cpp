#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
    const Eigen::MatrixXi & F,
    Eigen::SparseMatrix<double>& A)
{
    int V_size = F.maxCoeff() + 1;
    std::vector<std::vector<int>> boundaries;
    igl::boundary_loop(F, boundaries);

    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    int size = 0;

    for (int i = 0; i < boundaries.size(); i++) 
    {
        size += (boundaries[i]).size();
    }

    triplets.reserve(4 * size);

    for (int i = 0; i < boundaries.size(); i++) {
        std::vector<int> b = boundaries[i];

        for (int j = 0; j < b.size(); j++)
        {
            int v1 = b[j];
            int v2 = b[(j + 1) % b.size()];

            triplets.push_back(T(V_size + v2, v1, 0.5));
            triplets.push_back(T(v2, V_size + v1, -0.5));

            triplets.push_back(T(v1, V_size + v2, 0.5));
            triplets.push_back(T(V_size + v1, v2, -0.5));
        }
    }

    A.resize(V_size * 2, V_size * 2);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

