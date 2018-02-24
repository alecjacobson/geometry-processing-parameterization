#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
    std::vector<std::vector<int> >boundaries;
    igl::boundary_loop(F, boundaries);
    int counter = 0;
    for (int i = 0; i < boundaries.size();i ++) {
        counter += boundaries[i].size();
    }

    int numV = F.maxCoeff() + 1;
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(counter*4);
    
    int Val;
    for (int curBoundary =0 ; curBoundary < boundaries.size(); curBoundary ++) {
        int numB = boundaries[curBoundary].size();
        for (int i = 0; i < numB; i ++) {
            if (i+1 == numB) {
                Val = 0;
            }
            else {
                Val = i + 1;
            }
            tripletList.push_back(T(boundaries[curBoundary][i], numV + boundaries[curBoundary][Val], 0.5));
            tripletList.push_back(T(boundaries[curBoundary][Val], numV + boundaries[curBoundary][i], -0.5));
            
            tripletList.push_back(T(numV + boundaries[curBoundary][Val],boundaries[curBoundary][i], 0.5));
            tripletList.push_back(T(numV + boundaries[curBoundary][i],boundaries[curBoundary][Val],  -0.5));

        
        }
    }
    A.resize(2*numV, 2*numV);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    
}

