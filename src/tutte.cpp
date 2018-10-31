#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>

void tutte(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & U)
{
    Eigen::SparseMatrix<double> L, Aeq;
    Eigen::VectorXd Beq;
    Eigen::VectorXi boundary;
    Eigen::MatrixXd UV;
    igl::min_quad_with_fixed_data<double> data;

    igl::cotmatrix(V, F, L);
    igl::boundary_loop(F, boundary);
    igl::map_vertices_to_circle(V, boundary, UV);
    igl::min_quad_with_fixed_precompute(L, boundary, Aeq, false, data);
    
    // Linear term is zero
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(data.known.size() + data.unknown.size(), 2);
    igl::min_quad_with_fixed_solve(data, B, UV, Beq, U);
}

