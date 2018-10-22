#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
#include <igl/vector_area_matrix.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Replace with your code
  int V_size = F.maxCoeff()+1,L_size;
  A.resize(V_size*2,V_size*2);
  std::vector<std::vector<int>> L;
  igl::boundary_loop(F,L);
  typedef Eigen::Triplet<double> T;
  std::vector<T> list;
  for (int j=0; j<L.size(); j++){
  	L_size=L[j].size();
  for (int i=0; i<L_size; i++){
    list.push_back(T(L[j][i],V_size+L[j][(i+1)%L_size],0.25));
    list.push_back(T(V_size+L[j][(i+1)%L_size],L[j][i],0.25));
    list.push_back(T(V_size+L[j][i],L[j][(i+1)%L_size],-0.25));
    list.push_back(T(L[j][(i+1)%L_size],V_size+L[j][i],-0.25));
  }
  }
  A.resize(2*V_size,2*V_size);
  A.setFromTriplets(list.begin(), list.end()); 
  //igl::vector_area_matrix(F,A);
}

