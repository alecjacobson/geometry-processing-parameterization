#include "tutte.h"
#include <igl/cotmatrix.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/min_quad_with_fixed.h>

void tutte(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::MatrixXd & U)
{
	Eigen::VectorXi boundary;
	igl::boundary_loop(F, boundary);

	Eigen::MatrixXd UV;
	igl::map_vertices_to_circle(V, boundary, UV);

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	
	Eigen::SparseMatrix<double> Aeq;
	Eigen::MatrixXd B, Beq;
	B.resize(V.rows(), 2);
	igl::min_quad_with_fixed(L, B, boundary, UV, Aeq, Beq, false, U);
	
	// igl::min_quad_with_fixed_precompute(L, boundary, Aeq, false, data);
	// igl::min_quad_with_fixed_solve(data, B, UV, Beq, U);
}
