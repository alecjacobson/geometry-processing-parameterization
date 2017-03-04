#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <vector>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Get boundary vertices
  Eigen::VectorXi boundary_vertices; 
  igl::boundary_loop(F, boundary_vertices);
  int num_vertices = V.rows();
  int num_boundary = boundary_vertices.rows();
  
  // Map the boundary vertices to unit circle
  Eigen::MatrixXd circle(num_boundary, 2);
  igl::map_vertices_to_circle(V, boundary_vertices, circle);
  
  // Construct laplacian
  Eigen::SparseMatrix<double> L(num_vertices, num_vertices);
  igl::cotmatrix(V, F, L);

  // Setup inputs
  Eigen::VectorXd B(num_vertices), x(num_vertices), y(num_vertices);
  Eigen::SparseMatrix<double> Aeq;
  Eigen::MatrixXd Beq;
  igl::min_quad_with_fixed(L, B, boundary_vertices, circle.col(0), Aeq, Beq, false, x);
  igl::min_quad_with_fixed(L, B, boundary_vertices, circle.col(1), Aeq, Beq, false, y);
  
  U.resize(num_vertices, 2);
  for (int i = 0; i < num_vertices; i++) {
    U(i,0) = x[i];
    U(i,1) = y[i];
  }
}

