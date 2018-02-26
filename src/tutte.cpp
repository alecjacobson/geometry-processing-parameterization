#include "tutte.h"
#include <map>
#include <set>
#include "igl/boundary_loop.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/cotmatrix.h"

typedef Eigen::Triplet<double> Triplet;

/**
* For a given starting vertex, indicated by index on a triangle 
* with 3 vertices indexed by {0, 1, 2}; returns the other indices.
*
* @param index the index of the starting vertex. Range of [0, 2].
* @param otherID indicates the other vertex index to return. Range of {0, 1}.
*/
int getOtherIndex(int index, int otherID) {
	if (index == 0) {
		int otherIndices[2] = {1, 2};
		return otherIndices[otherID];
	} else if (index == 1) {
		int otherIndices[2] = {0, 2};
		return otherIndices[otherID];
	} else {
		int otherIndices[2] = {0, 1};
		return otherIndices[otherID];
	}
}

void makeLaplacian(const Eigen::MatrixXi & F,
		Eigen::VectorXi boundaryIndices, int vertexCount,
		Eigen::SparseMatrix<double>& laplacian) {

	laplacian.resize(vertexCount, vertexCount);

	std::vector<Triplet> triplets;
	int faceCount = F.rows();
	for (int i = 0; i < faceCount; i++) {
		for (int j = 0; j < 3; j++) {
			int i1 = F(i, getOtherIndex(j, 0));
			int i2 = F(i, getOtherIndex(j, 1));

			triplets.push_back(Triplet(i1, i2, -1));
			triplets.push_back(Triplet(i2, i1, -1));
			triplets.push_back(Triplet(i1, i1, 1));
			triplets.push_back(Triplet(i2, i2, 1));
		}
	}

	laplacian.setFromTriplets(triplets.begin(), triplets.end());
}

void tutte(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		Eigen::MatrixXd & U) {
	
	// Figure out what the boundary is.
	Eigen::VectorXi boundaryIndices;
	igl::boundary_loop(F, boundaryIndices);

	// Map the boundary vertices to the circle.
	Eigen::MatrixXd boundaryMapping(boundaryIndices.rows(), 2);
	igl::map_vertices_to_circle(V, boundaryIndices, boundaryMapping);

	// Construct our matrix, Q. But I think the min_quad_ function just
	// expects one quarter of it.
	Eigen::SparseMatrix<double> L;
	makeLaplacian(F, boundaryIndices, V.rows(), L);

	// Based on the equation in the tutorial,
	// http://libigl.github.io/libigl/tutorial/tutorial.html#chapter5:parametrization
	// many of the terms are set to 0 to match the minimization problem 
	// we actually have.
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(V.rows(), 2);
	Eigen::SparseMatrix<double> Aeq(1, V.rows());
	Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(1, 2);
	igl::min_quad_with_fixed(L, B, boundaryIndices, boundaryMapping, Aeq, Beq, true, U);
}

