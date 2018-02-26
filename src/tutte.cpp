#include "tutte.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  int vetr_num = V.rows();

  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::VectorXi BL;
  igl::boundary_loop(F, BL);

  Eigen::MatrixXd UV;
  igl::map_vertices_to_circle(V, BL, UV);

  bool pd = false;
  Eigen::VectorXd B = Eigen::VectorXd::Zero(vetr_num);
  Eigen::SparseMatrix<double> Aeq;
  Eigen::VectorXd Beq;
  igl::min_quad_with_fixed(L, B, BL, UV, Aeq, Beq, pd, U);
}
