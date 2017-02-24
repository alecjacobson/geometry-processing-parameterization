#include "tutte.h"

#include <vector>

#include "igl/boundary_loop.h"
#include "igl/cotmatrix.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/min_quad_with_fixed.h"

void tutte(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::MatrixXd & U)
{
  // Replace with your code
	U = V.leftCols(2);

	//First, map boundary to unit circle, and fix those points.
	Eigen::VectorXi boundary;
	igl::boundary_loop(F, boundary);
	Eigen::MatrixXd mapped_points;
	igl::map_vertices_to_circle(V, boundary, mapped_points);

	//Then, with those points fixed, construct the rest of the minimization problem.
	Eigen::SparseMatrix<double> cot;
	igl::cotmatrix(V, F, cot);

	//There are two optimizations to solve, really - one for each column of U. 
	Eigen::VectorXd Us = V.col(0);
	Eigen::VectorXd Vs = V.col(1);

	int n = cot.rows();
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, 1);
	Eigen::SparseMatrix<double> Aeq = Eigen::SparseMatrix<double>(0, 0);
	Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(0, 0);

	igl::min_quad_with_fixed_data<double> mqwf;
	igl::min_quad_with_fixed_precompute(cot, boundary, Aeq, false, mqwf);
	igl::min_quad_with_fixed_solve(mqwf, B, mapped_points.col(0), Beq, Us);
	igl::min_quad_with_fixed_solve(mqwf, B, mapped_points.col(1), Beq, Vs);	

	U.col(0) = Us;
	U.col(1) = Vs;
}

