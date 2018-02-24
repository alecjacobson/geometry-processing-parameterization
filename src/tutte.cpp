#include "tutte.h"
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  Eigen::VectorXi bnd;
  igl::boundary_loop(F, bnd);

  Eigen::MatrixXd bnd_mapping;
  igl::map_vertices_to_circle(V, bnd, bnd_mapping);

  Eigen::SparseMatrix<double> laplacian(V.rows(), V.rows());
  igl::cotmatrix(V, F, laplacian);

  Eigen::SparseMatrix<double> Aeq;
  Eigen::VectorXd Beq;
  Eigen::VectorXd B(V.rows(),1);
  B.setZero();

  U.resize(V.rows(),2);
  Eigen::VectorXd U_1 = U.col(0); 
  Eigen::VectorXd U_2 = U.col(1);

  igl::min_quad_with_fixed_data <double>mqwf;
  igl::min_quad_with_fixed_precompute(laplacian, bnd, Aeq, false, mqwf);
  igl::min_quad_with_fixed_solve(mqwf, B, bnd_mapping, Beq, U);

}

