#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <Eigen/SparseCholesky>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int n = V.rows();
  
  Eigen::SparseMatrix<double> L(n, n);
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  for (int f = 0; f < F.rows(); f++) {
    for (int e = 0; e < 3; e++) {
      int i = F(f, e);
      int j = F(f, (e+1)%3);
      double wij = 1 / (V.row(i) - V.row(j)).norm();
      tripletList.push_back(T(i, j, -wij));
      tripletList.push_back(T(j, i, -wij));
      tripletList.push_back(T(i, i, wij));
      tripletList.push_back(T(j, j, wij));
    }
  }
  
  L.setFromTriplets(tripletList.begin(), tripletList.end());
  
  Eigen::VectorXi bd;
  igl::boundary_loop(F, bd); // bd = longest boundary loop
  bd = bd.rowwise().reverse();

  Eigen::MatrixXd C;
  igl::map_vertices_to_circle(V, bd, C);
  
  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(L, bd, Eigen::SparseMatrix<double>(0, 0), false, data);
  
  Eigen::MatrixXd U1, U2;
  igl::min_quad_with_fixed_solve(data, Eigen::MatrixXd::Zero(n, 1), C.col(0), Eigen::MatrixXd::Zero(0,0), U1);
  igl::min_quad_with_fixed_solve(data, Eigen::MatrixXd::Zero(n, 1), C.col(1), Eigen::MatrixXd::Zero(0,0), U2);
  
  U.resize(n, 2);
  U << U1, U2;
  
  // for some reason U always has the wrong orientation, so we flip it along the x-axis...
  Eigen::Matrix2d R;
  R << 1, 0, 0, -1;
  U = U * R.transpose();
  // ...and pretend everything's OK
}

