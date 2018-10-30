#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

typedef Eigen::Triplet<double> T;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  int V_size = F.maxCoeff()+1;
  A.resize(V_size*2,V_size*2);
  std::vector<std::vector<int>> B;
  igl::boundary_loop(F,B);
  std::vector<T> tlist;

  for(int i = 0; i < B.size(); i++)
  {
    for(int j = 0; j < B[i].size(); j++){
    //Get adjacent
    int k = j+1;
    //Prevent out of range
    if (j == B[i].size() - 1){
      k = 0;
    }
    tlist.push_back(T(B[i][j], B[i][k] + V_size, 0.5));
    tlist.push_back(T(B[i][k], B[i][j] + V_size, -0.5));
    }
  }

  //Symmetric
  A.setFromTriplets(tlist.begin(), tlist.end());
  Eigen::SparseMatrix<double> AT;
  AT = A.transpose();
  A = A + AT;
  A = A * 0.5;
}

