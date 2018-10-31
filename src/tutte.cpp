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
	// igl::cotmatrix(V, F, L);
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	for (int f = 0; f < F.rows(); f++) {
		for (int i = 0; i < 3; i++) {
			int j = (i + 1) % 3;
			int k = (i + 2) % 3;
			double norm_ij = (V.row(F(f, i)) - V.row(F(f, j))).norm();
			double norm_ik = (V.row(F(f, i)) - V.row(F(f, k))).norm();

			// case i != j
			tripletList.push_back(T(F(f, i), F(f, j), 1.0/norm_ij));
			tripletList.push_back(T(F(f, i), F(f, k), 1.0/norm_ik));
			// case i == j
			tripletList.push_back(T(F(f, i), F(f, i), -1.0/norm_ij));
			tripletList.push_back(T(F(f, i), F(f, i), -1.0/norm_ik));
		}
	}
	L.resize(V.rows(), V.rows());
	L.setFromTriplets(tripletList.begin(), tripletList.end());
	
	Eigen::SparseMatrix<double> Aeq;
	Eigen::VectorXd B = Eigen::VectorXd::Zero(V.rows());
	igl::min_quad_with_fixed(L, B, boundary, UV, Aeq, B, false, U);
	
	// igl::min_quad_with_fixed_precompute(L, boundary, Aeq, false, data);
	// igl::min_quad_with_fixed_solve(data, B, UV, Beq, U);
}
