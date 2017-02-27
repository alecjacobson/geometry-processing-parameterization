#include "tutte.h"
#include "igl/boundary_loop.h"
#include "igl/cotmatrix.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/massmatrix.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/repdiag.h"
#include "igl/eigs.h"

using namespace Eigen;
using namespace igl;

typedef Eigen::Triplet<double> T;

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	VectorXi bVert;
	boundary_loop(F, bVert);

	MatrixXd UV;
	map_vertices_to_circle(V, bVert, UV);

	SparseMatrix<double> W;
	cotmatrix(V, F, W);
		
	SparseMatrix<double> Aeq(0,0); //dummy data
	min_quad_with_fixed_data<double> data;
	bool result = min_quad_with_fixed_precompute(W, bVert, Aeq, false, data);
	assert(result);

	int n = V.rows();
	VectorXd B(n);
	B.setZero();

	VectorXd Beq(0); //dummy data
	Eigen::VectorXd Z1(n);
	Eigen::VectorXd Z2(n);
	min_quad_with_fixed_solve(data, B, UV.col(0), Beq, Z1);
	min_quad_with_fixed_solve(data, B, UV.col(1), Beq, Z2);
	
	U.resize(n, 2);
	U.col(0) = Z1;
	U.col(1) = Z2;
}

