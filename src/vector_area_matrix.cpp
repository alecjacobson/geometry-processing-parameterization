#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  // Get boundary vertices
  std::vector<std::vector<int>> boundary_vertices; 
  igl::boundary_loop(F, boundary_vertices);
    
  // Construct A
  int vertex_num = F.maxCoeff() + 1;
  A.resize(2 * vertex_num, 2 * vertex_num);
  std::vector<Eigen::Triplet<double>> tripletList;
  
  std::vector<int> boundary;
  int vi, vj;
  for (int c = 0; c < boundary_vertices.size(); c++) {
    boundary = boundary_vertices[c];
    for (int i = 0; i < boundary.size() - 1; i++) {
      vi = boundary[i];
      vj = boundary[i+1];
      tripletList.push_back(Eigen::Triplet<double>(vertex_num + vi, vj, -0.5));
      tripletList.push_back(Eigen::Triplet<double>(vi, vertex_num + vj, 0.5));
      tripletList.push_back(Eigen::Triplet<double>(vj, vertex_num + vi, -0.5));
      tripletList.push_back(Eigen::Triplet<double>(vertex_num + vj, vi, 0.5));
    }
    
    // add another tripletList for the edge connecting the last vertex with the first
    vi = boundary[boundary.size() - 1];
    vj = boundary[0];
    tripletList.push_back(Eigen::Triplet<double>(vertex_num + vi, vj, -0.5));
    tripletList.push_back(Eigen::Triplet<double>(vi, vertex_num + vj, 0.5));
    tripletList.push_back(Eigen::Triplet<double>(vj, vertex_num + vi, -0.5));
    tripletList.push_back(Eigen::Triplet<double>(vertex_num + vj, vi, 0.5));  
  }

  A.setFromTriplets(tripletList.begin(), tripletList.end());
}