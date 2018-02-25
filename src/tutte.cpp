#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // // Replace with your code
  // U = V.leftCols(2);

	// Initialize the U container
	// U.resize(V.rows(), 2);

	// Strategy:
	// First, find all the edges of the mesh that lie on a boundary. We can use igl for this.
	// Second, use another whitelisted igl function to map the boundary edges to a circle.
	// Third, use the above as constraints to optimize for the remaining entries of U. Use igl to do this optimization.

	// Obtain boundary vertices
	std::vector<std::vector<int>> V_bdry_all;
	igl::boundary_loop(F, V_bdry_all);

	// We care about the outermost boundary, which is likely the longest one. So find the longest boundary.
	int len = 0;
	std::vector<int> V_bdry_vec;

	for (int ii = 0; ii < V_bdry_all.size(); ii++)
	{
		if (V_bdry_all[ii].size() > len)
		{
			len = V_bdry_all[ii].size();
			V_bdry_vec = V_bdry_all[ii];
		}
	}

	// Convert to Eigen vector
	Eigen::VectorXi V_bdry (V_bdry_vec.size());
	for (int ii = 0; ii < V_bdry_vec.size(); ii++)
		V_bdry(ii) = V_bdry_vec[ii];

	// I'm going to assume there's only one boundary

	// Now map these on to a circle
	Eigen::MatrixXd U_bdry;
	igl::map_vertices_to_circle(V, V_bdry, U_bdry);

	// Create the system that we're trying to minimize, with the cotangent matrix-based modification
	// We're going to use the igl function that minimizes Z'AZ, where A is the cotangent matrix.

	Eigen::SparseMatrix<double> A;
	igl::cotmatrix(V, F, A);

	// Create empty matrices for the equality constraints, of which we have none (because the known vertices are separate inputs)
	Eigen::SparseMatrix<double> Aeq;

	// Output container
	igl::min_quad_with_fixed_data<double> data;

	// Factorize
	igl::min_quad_with_fixed_precompute(A, V_bdry, Aeq, false, data);

	// Solve. The linear part of the cost function is just zero.
	Eigen::MatrixXd B (A.rows(), 2), Beq;
	B.setZero(B.rows(), 2);

	igl::min_quad_with_fixed_solve(data, B, U_bdry, Beq, U);

	return;

}

