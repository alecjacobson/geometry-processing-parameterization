#include "tutte.h"

#include <vector>
#include <igl/edges.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>

#include <iostream>
using std::cout;

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    // Fill L

    int n = V.rows();
    Eigen::SparseMatrix<double> L;
    std::vector<Eigen::Triplet<double>> triplets;

    double w;
    int i, j, k;
    for (int a = 0; a < F.rows(); ++a) {
        for (int b = 0; b < 3; ++b) {
            i = F(a, b);
            j = F(a, (b+1)%3);
            w = 1/(V.row(i) - V.row(j)).squaredNorm();
            triplets.emplace_back(i, j,  w);
            triplets.emplace_back(j, i,  w);
            triplets.emplace_back(i, i, -w);
            triplets.emplace_back(j, j, -w);
        }
    }

    L.resize(n, n);
    L.setFromTriplets(triplets.begin(), triplets.end());

    // Find 2D coordinate `UV` for boundary edges fixed to unit circle

    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);
    Eigen::MatrixXd UV;
    igl::map_vertices_to_circle(V, bnd, UV);

    // Minimize trace of a quadratic

    Eigen::SparseMatrix<double> Aeq;
    Eigen::MatrixXd Beq;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, 2);
    igl::min_quad_with_fixed_data<double> tutte_data;
    igl::min_quad_with_fixed_precompute(L, bnd, Aeq, false, tutte_data);
    igl::min_quad_with_fixed_solve(tutte_data, B, UV, Beq, U);
}

