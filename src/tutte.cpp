#include "tutte.h"

#include "igl/boundary_loop.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/cotmatrix.h"
#include "igl/min_quad_with_fixed.h"

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  Eigen::VectorXi boundary_verts;
  igl::boundary_loop(F, boundary_verts);

  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::MatrixXd boundary_pos_on_circ;
  igl::map_vertices_to_circle(V, boundary_verts, boundary_pos_on_circ);

  Eigen::SparseMatrix<double> sparse_zero;
 	Eigen::VectorXd zero = Eigen::VectorXd::Zero(V.rows());
  igl::min_quad_with_fixed(L, zero, boundary_verts, boundary_pos_on_circ, sparse_zero, zero, false, U);
}
