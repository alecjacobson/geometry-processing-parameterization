#include "vector_area_matrix.h"

#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
    // Convert asymetric Ã -> symmetric A
    // A = 1/2 (Ã + Ã')
    //  where Ã is asymetric area matrix
    //  s.t. for each (i,j) \in bnd(S),
    //      Ã[i, j+n] = 1
    //      Ã[i+n, j] = 1

    std::vector<std::vector<int>> bnds;
    igl::boundary_loop(F, bnds);

    int n = F.maxCoeff() + 1;
    std::vector<Eigen::Triplet<double>> triplets;
    int i, j;

    for (int a = 0; a < bnds.size(); ++a) {
        auto bnd = bnds[a];
        for (int b = 0; b < bnd.size(); ++b) {
            i = bnd[b];
            j = bnd[(b+1)%bnd.size()];
            triplets.emplace_back(i, j+n,  0.5);
            triplets.emplace_back(i+n, j, -0.5);
            triplets.emplace_back(j+n, i,  0.5);
            triplets.emplace_back(j, i+n, -0.5);
        }
    }

    A.resize(2*n, 2*n);
    A.setFromTriplets(triplets.begin(), triplets.end());
}