#include "vector_area_matrix.h"
#include <iostream>
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // // Replace with your code
  // A.resize(F.maxCoeff()+1,F.maxCoeff()+1);

	// First figure out number of vertices
	int N_vert = F.maxCoeff()+1;

	// This matrix should be a bunch of 1s and -1s corresponding to boundary vertices.

	// Obtain boundary vertices
	std::vector<std::vector<int>> V_bdry_all;
	igl::boundary_loop(F, V_bdry_all);

	// Initialize the containers to store triplets, to assemble the sparse area matrix. This is based on Eigen documentation.
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets_for_A, triplets_for_AT;
	//tripletList.reserve(10*V_bdry.cols());
	
	// Now loop through all the boundaries that exist
	std::vector<int> V_bdry;
	for (int ii = 0; ii < V_bdry_all.size(); ii++)
	{
		// Extract this boundary
		V_bdry = V_bdry_all[ii];

		// Now loop through each vertex in this boundary
		for (int jj = 0; jj < V_bdry.size(); jj++)
		{
			int temp = jj + 1;
			if (temp == V_bdry.size())
				temp = 0;
			int index1 = V_bdry[(jj + 1) % V_bdry.size()];
			int index2 = V_bdry[jj];

			// Collect the values and indices for the area matrix A
			triplets_for_A.push_back(T(index2, index1 + N_vert, 1.0));
			triplets_for_A.push_back(T(index2 + N_vert, index1, -1.0));

			// Collect the values and indices for the area matrix transpose A^T by swapping indices
			triplets_for_AT.push_back(T(index1 + N_vert, index2, 1));
			triplets_for_AT.push_back(T(index1, index2 + N_vert, -1.0));
		}
	}

	// Assemble the matrices A and A^T
	Eigen::SparseMatrix<double> AT(2*N_vert, 2*N_vert);
	A.resize(2*N_vert, 2*N_vert);
	A.setFromTriplets(triplets_for_A.begin(), triplets_for_A.end());
	AT.setFromTriplets(triplets_for_AT.begin(), triplets_for_AT.end());

	// Symmetrize A and divide by 2
	A = 0.5*(A + AT);

	return;

}

