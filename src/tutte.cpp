#include "tutte.h"
#include <igl/cotmatrix.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    Eigen::SparseMatrix<double> Aeq, L(V.rows(), V.rows());
    Eigen::VectorXd B, Beq = Eigen::VectorXd::Zero(V.rows());
    Eigen::VectorXi bl;
    Eigen::MatrixXd UV;
    
    igl::cotmatrix(V, F, L);
    igl::boundary_loop(F, bl);
    igl::map_vertices_to_circle(V, bl, UV);
    igl::min_quad_with_fixed_data<double> data;
    igl::min_quad_with_fixed_precompute(L, bl, Aeq, false, data);
    igl::min_quad_with_fixed_solve<double>(data, B, UV, Beq, U);
    U.col(0) = -U.col(0);
}
