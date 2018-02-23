#include "vector_area_matrix.h"
#include <igl/boundary_facets.h>
#include <vector>
typedef Eigen::Triplet<double> tuple;
void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  //I don't fully understand the definition of matrix A in the README file,
  //so I use the official implementation of the igl::vector_area_matrix() as a template
  //During the testing, I discovered a bug in the official implementation and fixed it
  const int nV = F.maxCoeff() + 1;
  Eigen::MatrixXi E;
  igl::boundary_facets(F, E);
  std::vector<tuple> tuple_list;
  tuple_list.reserve(2*E.rows());
  for (int k = 0; k < E.rows(); k++)
  {
    tuple_list.push_back(tuple(E(k, 0) + nV, E(k, 1), 0.5));
    tuple_list.push_back(tuple(E(k, 1) + nV, E(k, 0), -0.5));
  }
  A.setFromTriplets(tuple_list.begin(), tuple_list.end());
  A += (Eigen::SparseMatrix<double>)A.transpose();
}

