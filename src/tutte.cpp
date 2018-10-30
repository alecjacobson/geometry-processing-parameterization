#include "tutte.h"
#include <iostream>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

using namespace std;

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  U = V.leftCols(2);
  // Detect the boundary vertices
  Eigen::VectorXi bnd;
  igl::boundary_loop(F,bnd);
  std::vector<std::vector<int>> lb;
  igl::boundary_loop(F, lb);

  // Map the boundary vertices to a circle
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V,bnd,bnd_uv);
  // Changing the orientation of the boudary circle
  bnd_uv.col(0) = -bnd_uv.col(0);

  // Calculate L
  Eigen::SparseMatrix<double> L(V.rows(), V.rows());
  L.setZero();
  std::vector<Eigen::Triplet<double>> triplets;
  // Adding weights for each half edge
  for (int f = 0; f < F.rows(); f++) {
      int i = F(f, 1);   int j = F(f, 2);
      double wij = 1.0/ (V.row(i) - V.row(j)).norm();
      //wij = 1;
      triplets.push_back({i, j,  wij });
      triplets.push_back({j, i,  wij });
      triplets.push_back({i, i, -wij });
      triplets.push_back({j, j, -wij });

      i = F(f, 2);   j = F(f, 0);
      triplets.push_back({i, j,  wij });
      triplets.push_back({j, i,  wij });
      triplets.push_back({i, i, -wij });
      triplets.push_back({j, j, -wij });

      i = F(f, 0);   j = F(f, 1);
      triplets.push_back({i, j,  wij });
      triplets.push_back({j, i,  wij });
      triplets.push_back({i, i, -wij });
      triplets.push_back({j, j, -wij });
  }
  L.setFromTriplets(triplets.begin(), triplets.end());

  // Calculate minimizer subject to the fixed boundary points
  igl::min_quad_with_fixed_data<double> data;
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(L, bnd, Aeq, false, data);

  U = Eigen::MatrixXd::Zero(data.n,2);
  Eigen::MatrixXd B(data.n, data.n);
  B.setZero();
  Eigen::MatrixXd Beq;
  igl::min_quad_with_fixed_solve(data, B, bnd_uv, Beq, U);
  igl::min_quad_with_fixed(L, B, bnd, bnd_uv, Aeq, Beq, false, U);
}
