#include "vector_area_matrix.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
    int n = F.maxCoeff()+1;
    
    std::vector<std::vector<int> > Loop;
    igl::boundary_loop(F,Loop);
    
    A.resize(2*n,2*n);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(int i = 0; i < Loop.size(); i++){
        for(int j = 0; j < Loop.at(i).size(); j++){
            int endpointa = Loop.at(i).at(j);
            int endpointb;
            if(j+1 == Loop.at(i).size()){
                endpointb = Loop.at(i).at(0);
            } else {
                endpointb = Loop.at(i).at(j+1);
            }
            tripletList.push_back(T(endpointa,endpointb + n,0.5));
            tripletList.push_back(T(endpointa + n,endpointb,-0.5));
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A = 0.5*(A + (Eigen::SparseMatrix<double>)A.transpose());
}

