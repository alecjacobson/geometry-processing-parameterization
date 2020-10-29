#include "tutte.h"
#include "igl/boundary_loop.h"
#include "igl/repdiag.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/cotmatrix.h"
#include "igl/min_quad_with_fixed.h"

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    Eigen::VectorXi b;
    igl::boundary_loop(F,b);

    Eigen::MatrixXd bc;
    igl::map_vertices_to_circle(V,b,bc);

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);

    Eigen::SparseMatrix<double> Q = -L;

    igl::min_quad_with_fixed_data<double> data;
    igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), false, data);

    const int dims = b.cols();
    const int size = data.n;

    igl::min_quad_with_fixed_solve(data, Eigen::MatrixXd::Zero(size, dims),
                                     bc, Eigen::MatrixXd::Zero(0, dims), U);
}

