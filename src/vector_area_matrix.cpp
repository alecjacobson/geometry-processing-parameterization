#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A) {
    
    // get the list of boundary loops
    std::vector<std::vector<int>> bound;
    igl::boundary_loop(F, bound);
    
    // build a triplet list to populate A.
    typedef Eigen::Triplet<double> T;
    std::vector<T> list;
    
    // cycle through each boundary loop to populate A
    int n = F.maxCoeff() + 1;
    for (int i = 0; i < bound.size(); i++) {
        
        std::vector<int> b = bound[i];
        
        for (int j = 0; j < b.size(); j++) {
            
            // get the next pair of adjacent vertices in the loop
            // (modulo is used so that the last vertex will circle back
            // to the first one.)
            int v0 = b[j];
            int v1 = b[(j + 1) % b.size()];
        
            // populate the coefficients matrix: because the naive construction
            // is not necessarily symmetric, we symmetrize our matrix by summing
            // the naive matrix with its transpose and multiplying by 0.5.
            
            // insert the elements from A
            list.push_back(T(v0 + n, v1, -0.5));
            list.push_back(T(v1 + n, v0, 0.5));
        
            // add the elements from A transpose
            list.push_back(T(v1, v0 + n, -0.5));
            list.push_back(T(v0, v1 + n, 0.5));
        }
    }
    
    A.resize(2*n, 2*n);
    A.setFromTriplets(list.begin(), list.end());
}