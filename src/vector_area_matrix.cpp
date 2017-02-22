#include "vector_area_matrix.h"
#include <vector>
#include <iostream>


// [a c] A [b d]^* = ad - bc
// [a c] [0 1; 2 3] [b d]^* = ad - bc
// [a c] [0b 1d; 2b 3d] = ad - bc
// [a0b a1d; c2b c3d] = ad - bc
// [_ +; - _] = ad - bc
// 0,1=1, 1,0=-1



void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{


  int dim = F.maxCoeff()+1;
  A = Eigen::SparseMatrix<double>(2*dim,2*dim);

  std::vector<Eigen::Triplet<double>> trips;


  for(int t = 0; t < F.rows(); ++t) {
      const auto tri = F.row(t);
      for(int i = 0; i < 3; ++i) {
          int a = i;
          int b = (i+1)%3;

          a = tri(a);
          b = tri(b);
          trips.emplace_back(a,b+dim,1);
          trips.emplace_back(a+dim,b,-1);

      }
  }

  A.setFromTriplets(trips.begin(),trips.end());


}

