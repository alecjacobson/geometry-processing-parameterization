#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

typedef Eigen::Triplet<double> Triplet;

void vector_area_matrix(
		const Eigen::MatrixXi & F,
		Eigen::SparseMatrix<double>& A) {
	
	// Our matrix depends soley on the boundary vertices.
	std::vector<std::vector<int>> boundaryLoops;
	igl::boundary_loop(F, boundaryLoops);

	// The dimensions of A are 2 * vertexCount to fit both axes of data.
	int vertexCount = F.maxCoeff() + 1;
	A.resize(2 * vertexCount, 2 * vertexCount);
	std::vector<Triplet> triplets;
	for (auto boundary : boundaryLoops) {
		for (int i = 0; i < boundary.size(); i++) {

			int indexA = boundary[i];
			int indexB = boundary[((i + 1) % boundary.size())];

			assert(indexA < vertexCount);
			assert(indexB < vertexCount);

			// This arrangment ensures that A represents the transformation
			// needed to compute determinants.
			triplets.push_back(Triplet(indexA, indexB + vertexCount, 0.5));
			triplets.push_back(Triplet(indexA + vertexCount, indexB, -0.5));
		}
	}
	A.setFromTriplets(triplets.begin(), triplets.end());

	// According to the ReadMe, what we had was A-bar. This forms A-actual.
	// Had to put in cast to avoid Eigen error.
	A = A + (Eigen::SparseMatrix<double>)A.transpose();
}

