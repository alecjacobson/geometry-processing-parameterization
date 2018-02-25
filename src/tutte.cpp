#include "tutte.h"
#include <igl/edges.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>

void uniform_laplacian(
	const Eigen::MatrixXi & E,
	Eigen::SparseMatrix<double> & L)
{
	// compute uniform laplacian for index i != j
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(E.rows() * 2);
	for (int i = 0; i < E.rows(); i++) {
		triplets.push_back(T(E(i, 0), E(i, 1), 1.0));
		triplets.push_back(T(E(i, 1), E(i, 0), 1.0));
	}
	int L_size = E.maxCoeff() + 1;
	L.resize(L_size, L_size);
	L.setFromTriplets(triplets.begin(), triplets.end());

	// compute diagonal part (−∑_{l≠i}L_{il}) where i = j
	for (int i = 0; i < L_size; i++) {
		L.insert(i, i) = -L.row(i).sum();
	}
}

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
	// Compute graph laplacian
	Eigen::MatrixXi E;
	igl::edges(F, E);
	Eigen::SparseMatrix<double> L;
	uniform_laplacian(E, L);

	// Compute ordered index of vertices of longest boundary loop
	Eigen::VectorXi bndV;
	igl::boundary_loop(F, bndV);

	// Map longest boundary vertices to unit circle
	Eigen::MatrixXd UV;
	igl::map_vertices_to_circle(V, bndV, UV);

	// Energy Minimization with boundary vertices mapped to unit circle constraint
	// Use min_quad_with_fixed with "bndV mapped to UV" constraint, Aeq*U = Beq is not used
	Eigen::SparseMatrix<double> Aeq;
	Eigen::VectorXd Beq;
	Eigen::VectorXd B = Eigen::VectorXd::Zero(L.rows());
	igl::min_quad_with_fixed(L, B, bndV, UV, Aeq, Beq, false, U); // Set to false since L is negative semi-definite not positive definite
}