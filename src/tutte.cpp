#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/cotmatrix.h>
#include <vector>
#include <iostream>

int get_index(const int num, const Eigen::VectorXi &arr) {
  int index = -1;
  for (int i = 0; i < arr.rows(); i++) {
    if (arr[i] == num) {
      index = i;
    }
  }
  return index;
}

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Get boundary vertices
  Eigen::VectorXi boundary_vertices; 
  igl::boundary_loop(F, boundary_vertices);
  int num_vertices = V.rows();
  int num_boundary = boundary_vertices.rows();
  
  // Get non-boundary vertices
  Eigen::VectorXi non_boundary_vertices(num_vertices - num_boundary);
  int counter = 0;
  for (int i = 0; i < num_vertices; i++) {
    if (get_index(i, boundary_vertices) == -1) {
      non_boundary_vertices[counter] = i;
      counter++;
    }
  }
  
  // Map the boundary vertices to unit circle
  Eigen::MatrixXd circle(num_boundary, 2);
  igl::map_vertices_to_circle(V, boundary_vertices, circle);
  
  // Construct laplacian
  Eigen::SparseMatrix<double> L(num_vertices, num_vertices);
  igl::cotmatrix(V, F, L);

  // We need to construct 2 subsets of the Laplacian
  Eigen::SparseMatrix<double> XX(num_vertices - num_boundary, num_vertices - num_boundary);
  Eigen::SparseMatrix<double> XY(num_vertices - num_boundary, num_boundary);
  std::vector<Eigen::Triplet<double>> tripletXX;
  std::vector<Eigen::Triplet<double>> tripletXY;
  
  int pb, qb, pn, qn;
  for (int i = 0; i < L.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L, i); it; ++it) {
      qb = get_index(it.col(), boundary_vertices);
      pn = get_index(it.row(), non_boundary_vertices);
      qn = get_index(it.col(), non_boundary_vertices);
      
      if ( pn > -1 && qb > -1 ) {
        tripletXY.push_back(Eigen::Triplet<double>(pn, qb, it.value()));
      } else if ( pn > -1 && qn > -1 ) {
        tripletXX.push_back(Eigen::Triplet<double>(pn, qn, it.value()));
      }
    }
  }
  XX.setFromTriplets(tripletXX.begin(), tripletXX.end());
  XY.setFromTriplets(tripletXY.begin(), tripletXY.end());

  // Now we can solve
  Eigen::VectorXd x, y;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
  cg.compute(XX);
  x = cg.solve(-1 * XY * circle.col(0));
  y = cg.solve(-1 * XY * circle.col(1));
  
  // Now we can build U
  U.resize(num_vertices, 2);
  for (int i = 0; i < num_boundary; i++) {
    U.row(boundary_vertices[i]) = circle.row(i);
  }
  for (int i = 0; i < num_vertices - num_boundary; i++) {
    U(non_boundary_vertices[i], 0) = x[i];
    U(non_boundary_vertices[i], 1) = y[i];
  }
}

