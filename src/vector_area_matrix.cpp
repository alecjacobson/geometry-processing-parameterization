#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  int V_size = F.maxCoeff()+1;
  A.resize(V_size*2,V_size*2);

  Eigen::VectorXi bnd;
  igl::boundary_loop(F,bnd);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tlist;
  for (int i=0; i<bnd.rows(); i++) {
    for (int j=0; j<bnd.cols(); j++) {
      if (j == bnd.cols() - 1) {
        tlist.push_back(T(bnd(i,j),bnd(i,0)+V_size,0.5));
        tlist.push_back(T(bnd(i,0), bnd(i,j)+V_size,-0.5));
        tlist.push_back(T(bnd(i,0)+V_size,bnd(i,j),0.5));
        tlist.push_back(T(bnd(i,j)+V_size,bnd(i,0),-0.5));
      }
      else {
        tlist.push_back(T(bnd(i,j),bnd(i,j+1)+V_size,0.5));
        tlist.push_back(T(bnd(i,j+1), bnd(i,j)+V_size,-0.5));
        tlist.push_back(T(bnd(i,j+1)+V_size,bnd(i,j),0.5));
        tlist.push_back(T(bnd(i,j)+V_size, bnd(i,j+1),-0.5));
      }
    }
  }
  A.setFromTriplets(tlist.begin(),tlist.end());
}

