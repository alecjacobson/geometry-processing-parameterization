#include "tutte.h"
#include <Eigen/Sparse>
#include <igl/edges.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>

using namespace std;

void graph_Laplacian(
  const Eigen::MatrixXi & E,
  Eigen::SparseMatrix<double> & L)
{
  typedef Eigen::Triplet<double> T;

  std::vector<T> tripletList;
  // For each edge, two elements of L are filled in with + 1
  // We will add the diagonal elements after
  tripletList.reserve(E.rows() * 2);

  for(int edge_number = 0; edge_number < E.rows(); edge_number++)
  {
    auto start_node_index = E(edge_number, 0);
    auto end_node_index = E(edge_number, 1);
    
    tripletList.push_back(T(start_node_index, end_node_index, 1.0));
    tripletList.push_back(T(end_node_index, start_node_index, 1.0));
  }
  L.setFromTriplets(tripletList.begin(), tripletList.end());
  
  // Set up Laplacian equality: what leaves a node, enters it
  for (int diagIndex = 0; diagIndex < L.rows(); diagIndex++){
    L.coeffRef(diagIndex, diagIndex) = -1 * L.row(diagIndex).sum();
  }
}

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int number_of_vertices = V.rows();
  
  // Borrowed from the libigl tutorial "Quadratic energy minimization" section

  ////// Obtain Q and B ////// 
  
  // The Graph Laplacian is our nxn matrix "Q" that appears in the quadratic form
  // Obtain graph edges "E"
  Eigen::MatrixXi E;
  igl::edges(F, E);
  // Construct the Graph Laplacian from edges
  Eigen::SparseMatrix<double> Q(number_of_vertices, number_of_vertices);
  graph_Laplacian(E, Q);
  // A Graph Laplacian always has at least one 0 eigen value, so it's not PD
  auto positive_definite = false;        
  // "B" is a n×1 vector of linear coefficients. We have none (all zero).
  Eigen::VectorXd B = Eigen::VectorXd::Zero(number_of_vertices);

  ////// Obtain Constraints ////// 

  // Obtain the loops involved in boundary conditions. 
  // In this case, we want to fix disk boundaries obtained from our mesh
  Eigen::VectorXi boundaries;
  igl::boundary_loop(F, boundaries);

  // "bc" corresponds to the values of the fixed vertices
  // In this case, we wish to map to the unit circle
  Eigen::MatrixXd bc;
  igl::map_vertices_to_circle(V, boundaries, bc);

  // Aeq is a is a sparse matrix of linear equality constraint coefficients 
  // Beq is a m×1 vector of linear equality constraint right-hand side values.
  // We have no such constraints.
  Eigen::SparseMatrix<double> Aeq;
  Eigen::VectorXd Beq;

  ////// Precomputation ////// 

  // mqwf contains the system matrix factorization 
  // it is used during solving with arbitrary linear terms, 
  // known values, and constraint in the right-hand sides
  igl::min_quad_with_fixed_data<double> mqwf;
  igl::min_quad_with_fixed_precompute(Q, boundaries, Aeq, positive_definite, mqwf);

  ////// Solve ////// 

  // Minimize the squared energy along the Laplacian
  igl::min_quad_with_fixed_solve<double>(mqwf, B, bc, Beq, U);
}

