#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
    
    //Define an array of arrays to store boundary vertices
    std::vector<std::vector<int> >boundaries;
    igl::boundary_loop(F, boundaries);
    
    //Count the number of boundary terms, this allows us to preallocate the right amount of space
    int counter = 0;
    for (int i = 0; i < boundaries.size();i ++) {
        counter += boundaries[i].size();
    }

    int numV = F.maxCoeff() + 1;
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(counter*4);
    
    int Val;
    
    //Loop through each boundary
    for (int curBoundary =0 ; curBoundary < boundaries.size(); curBoundary ++) {
        int numB = boundaries[curBoundary].size();
        
        //For each boundary, add a triplet for the cells of A and the cells of A^t
        for (int i = 0; i < numB; i ++) {
            if (i+1 == numB) {
                Val = 0;
            }
            else {
                Val = i + 1;
            }
            
            //It looks to me like igl::vector_area_matrix uses 0.25 here instead of 0.5.
            //This would means that [L 0; 0 L] should be halfed in the generalized eigensolver
            tripletList.push_back(T(boundaries[curBoundary][i], numV + boundaries[curBoundary][Val], 0.5));
            tripletList.push_back(T(boundaries[curBoundary][Val], numV + boundaries[curBoundary][i], -0.5));
            
            tripletList.push_back(T(numV + boundaries[curBoundary][Val],boundaries[curBoundary][i], 0.5));
            tripletList.push_back(T(numV + boundaries[curBoundary][i],boundaries[curBoundary][Val],  -0.5));

        
        }
    }
    A.resize(2*numV, 2*numV);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    
}

