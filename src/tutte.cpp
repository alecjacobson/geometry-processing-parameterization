#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/cotmatrix.h>
#include <vector>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
    //Get the boundary
    Eigen::VectorXi Lb;
    igl::boundary_loop(F, Lb);
    Eigen::MatrixXd boundaryUV; 
    igl::map_vertices_to_circle(V, Lb ,boundaryUV);

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    igl::min_quad_with_fixed_data<double> data;
    Eigen::SparseMatrix<double> dummy(0,0);
    igl::min_quad_with_fixed_precompute(L, Lb, dummy, false, data);
    Eigen::MatrixXd Z;
    Eigen::MatrixXd dummy2(0,0); 
    igl::min_quad_with_fixed_solve(data, Eigen::MatrixXd::Zero(V.rows(), 2),boundaryUV, dummy2, Z);
    U = Z;
}

