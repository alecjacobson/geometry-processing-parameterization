#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <vector>
#include <igl/cotmatrix.h>
void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Map boundary vertices to unit circle
  U = Eigen::MatrixXd::Zero(V.rows(), 2);
  Eigen::VectorXi bnd;
  igl::boundary_loop(F,bnd); // vertices of boundary
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V,bnd,bnd_uv);

  // Get Laplacian
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  // Linear term is 0 except at known values
  Eigen::VectorXd B = Eigen::VectorXd::Zero(V.rows(),1);

  // Empty constraints
  Eigen::VectorXd Beq;
  Eigen::SparseMatrix<double> Aeq;
  // Minimize energey trace( 0.5*U'*A*U + U'*B + constant )
  igl::min_quad_with_fixed_data<double> tutte;
  igl::min_quad_with_fixed_precompute(L, bnd, Aeq,false, tutte);
  igl::min_quad_with_fixed_solve(tutte, B, bnd_uv, Beq, U);
}
