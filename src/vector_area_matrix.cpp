#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
#include <iostream>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Get boundary vertices
  std::vector<int> boundary_vertices; 
  igl::boundary_loop(F, boundary_vertices);
  
  // Construct A
  int vertex_num = F.maxCoeff() + 1;
  A.resize(2*vertex_num, 2*vertex_num);
  A.setZero();
  std::vector<Eigen::Triplet<double>> tripletList;
  
  int vi, vj;
  for (int i = 0; i < boundary_vertices.size() - 1; i++) {
    vi = boundary_vertices[i];
    vj = boundary_vertices[i+1];
    std::cout << vi << std::endl;
    tripletList.push_back(Eigen::Triplet<double>(vi, vertex_num + vj, 1));
    tripletList.push_back(Eigen::Triplet<double>(vertex_num + vj, vi, 1));
    tripletList.push_back(Eigen::Triplet<double>(vertex_num + vi, vj, -1));
    tripletList.push_back(Eigen::Triplet<double>(vj, vertex_num + vi, -1));
  }
  
  // add another tripletList for the edge connecting the last vertex with the first
  vi = boundary_vertices[boundary_vertices.size() - 1];
  vj = boundary_vertices[0];
  tripletList.push_back(Eigen::Triplet<double>(vi, vertex_num + vj, 1));
  tripletList.push_back(Eigen::Triplet<double>(vertex_num + vi, vj, -1));
  
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  A = 0.5 * A;
}