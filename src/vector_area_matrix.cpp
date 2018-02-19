#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
using namespace std;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  // A.resize(F.maxCoeff()+1,F.maxCoeff()+1);

  // My code
	int numV = F.maxCoeff() + 1;

	// extract all boundary loops
	vector<vector<int> > allLoops;
	igl::boundary_loop(F,allLoops);

	// assemble triplet list
	Eigen::VectorXi bdLoop; 
	std::vector< Eigen::Triplet<double> > tripletList;
	tripletList.reserve(allLoops.size() * 5);
	for (int ii = 0; ii < allLoops.size(); ii++){
		// resize for this loop
		bdLoop.resize(allLoops[ii].size()); 
		
		// assign bd indices to bdloop
		for (int jj = 0; jj < allLoops[ii].size(); jj++){ 
			bdLoop(jj) = allLoops[ii][jj];
		}

		// add indices to triplet
		for (int ii = 0; ii < bdLoop.rows(); ii++){
			int VaIdx = bdLoop(ii);
			int VbIdx = bdLoop((ii+1) % bdLoop.rows());
			tripletList.emplace_back(VaIdx+numV, VbIdx, -0.5); 
			tripletList.emplace_back(VaIdx, VbIdx+numV,  0.5);
		}
	}

	// assemble A matrix
	A.resize(numV*2, numV*2);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A += Eigen::SparseMatrix<double>(A.transpose()); // force it to be symmetric
}

