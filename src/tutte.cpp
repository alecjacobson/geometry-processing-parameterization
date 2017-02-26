#include "tutte.h"
#include <igl/cotmatrix.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/min_quad_with_fixed.h>
#include <vector>


void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  U = V.leftCols(2);
	Eigen::VectorXi boundry;
	igl::boundary_loop(F, boundry);

	// Map boundry to circle
	Eigen::MatrixXd Y;
	igl::map_vertices_to_circle(V, boundry, Y);

	// Compute cotangent matrix
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	// Solve two linear systems. One for each parametric coordinate
	igl::min_quad_with_fixed_data<double> data;
	Eigen::SparseMatrix<double> A = -L;
	Eigen::SparseMatrix<double> Aeq;
	min_quad_with_fixed_precompute(A, boundry, Aeq, true, data);
	Eigen::VectorXd B = Eigen::VectorXd::Zero(V.rows(), 1);
	Eigen::VectorXd Beq;
	Eigen::VectorXd Y1 = Y.col(0);
	Eigen::VectorXd U1;
	Eigen::VectorXd Y2 = Y.col(1);
	Eigen::VectorXd U2;
	igl::min_quad_with_fixed_solve(data, B, Y1, Beq, U1);
	igl::min_quad_with_fixed_solve(data, B, Y2, Beq, U2);

	// Concatinate parametric coordinates
	U << U1, U2;
}

