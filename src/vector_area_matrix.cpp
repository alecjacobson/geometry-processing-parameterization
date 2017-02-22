#include "vector_area_matrix.h"
#include <vector>
#include <iostream>
#include <igl/boundary_loop.h>
#include <igl/vector_area_matrix.h>


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

  using Index = Eigen::MatrixXi::Scalar;
  std::vector<Index> bnd_loop;
  igl::boundary_loop(F,bnd_loop);

  for(int i = 0; i < bnd_loop.size(); ++i) {
      int a = bnd_loop[i];
      int b = bnd_loop[(i+1)%bnd_loop.size()];

      trips.emplace_back(a,b+dim,1);
      trips.emplace_back(a+dim,b,-1);
      trips.emplace_back(b+dim,a,1);
      trips.emplace_back(b,a+dim,-1);

  }

  A.setFromTriplets(trips.begin(),trips.end());
  A = .5 * A;


  igl::vector_area_matrix(F,A);
  //std::cout << "real: " << std::endl;
  //std::cout << 4*A << std::endl;
}

