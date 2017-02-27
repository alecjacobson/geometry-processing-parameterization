#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/cotmatrix.h>

#include <igl/harmonic.h>
void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	// Replace with your code
  
	//LibIGL solution to compare too...
	//Eigen::VectorXi bnd;
	//igl::boundary_loop(F, bnd);
	//Eigen::MatrixXd bnd_uv;
	//igl::map_vertices_to_circle(V, bnd, bnd_uv);
	//igl::harmonic(V, F, bnd, bnd_uv, 1, U);

	//Get the disk boundary of the mesh, and map it to a circle
	Eigen::VectorXi boundary;
	igl::boundary_loop(F, boundary);
	Eigen::MatrixXd mapping;
	igl::map_vertices_to_circle(V, boundary, mapping);

	//Find the cotangent matrix for the discrete minimization problem
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	//Solve the two systems corresponding to each coordinate U1 and U2
	
	//Set the linear equality constraints and linear coefficients 	
	Eigen::SparseMatrix<double> Aeq;
	Eigen::VectorXd B = Eigen::VectorXd::Zero(L.rows(), 1); // Define the Size of B or else get error
	Eigen::VectorXd Beq, U1, U2;

	igl::min_quad_with_fixed_data<double> data;
	igl::min_quad_with_fixed_precompute(L, boundary, Aeq, false, data);
	igl::min_quad_with_fixed_solve(data, B, mapping.col(0), Beq, U1);
	igl::min_quad_with_fixed_solve(data, B, mapping.col(1), Beq, U2);

	U = Eigen::MatrixXd(V.rows(), 2);
	U.col(0) = U1;
	U.col(1) = U2;




}

