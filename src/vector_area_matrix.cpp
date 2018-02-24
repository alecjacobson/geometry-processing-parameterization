#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"
#include <iostream>

using namespace std;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int num_nodes = F.maxCoeff() + 1;

  // When we solve, U [n, 2] is flattened --> U = vector of dimension n*2
  A.resize(2 * num_nodes, 2 * num_nodes);

  // Obtain boundaries of mesh, each of which are closed loops
  std::vector<std::vector<int>> Mesh_Loops;
  igl::boundary_loop(F, Mesh_Loops);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  // For each boundary
  // i = boundary number
  for(int i = 0; i < Mesh_Loops.size(); i++){
    // For each point in the boundary (consecutive and oriented in one direction)
    // j = point number
    for(int j = 0; j < Mesh_Loops[i].size(); j++){
        int start_vertex_index = Mesh_Loops[i][j];
        int end_vertex_index;

        if (j == Mesh_Loops[i].size() - 1)
        {
          // End point, loop back to start point
          end_vertex_index = Mesh_Loops[i][0];
        } 
        else 
        {
          end_vertex_index = Mesh_Loops[i][j+1];
        }

        tripletList.push_back(T(start_vertex_index, end_vertex_index + num_nodes, 1.0));
        tripletList.push_back(T(start_vertex_index + num_nodes, end_vertex_index, -1.0));
    }
  }

  // Set elements
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  // Constant of a half occurs in the integration of interpolated points
  A = 0.5 * A;

  // Symmetrize A using common trick
  Eigen::SparseMatrix<double> A_T = A.transpose();
  A = 0.5 * (A + A_T);
}

