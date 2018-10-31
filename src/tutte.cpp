#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Known coefficients: Boundary indexes and mapped vertices
  Eigen::VectorXi bound_index;
  Eigen::MatrixXd bound_vertices;
  igl::boundary_loop(F, bound_index);
  igl::map_vertices_to_circle(V, bound_index, bound_vertices);

  // Quadratic coefficients matrix: Contagent Laplacian
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // Precompute quadratic functions:
  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(L, bound_index, 
                Eigen::SparseMatrix<double>(), false, data);
  
  // Solve quardratic functions:
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(data.n, 2);
  igl::min_quad_with_fixed_solve(data, B, bound_vertices, 
                Eigen::MatrixXd(), U);
}

