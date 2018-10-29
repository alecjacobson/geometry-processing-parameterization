#include "tutte.h"
#include "igl/boundary_loop.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/edges.h"
#include "igl/min_quad_with_fixed.h"
#include <iostream>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  
  // Determine all boundary loop on the mesh
  std::vector<std::vector<int>> bdL;
  igl::boundary_loop(F, bdL);

  // Identify the largest boundary loop
  int max_idx = 0;
  for (int i = 0; i < bdL.size(); i++) {
	if (bdL[i].size() >= bdL[max_idx].size()) {
		max_idx = i;
	}
  } 
  Eigen::VectorXi bdry(bdL[max_idx].size());
  for (int i = 0; i < bdL[max_idx].size(); i++) {
	bdry[i] = bdL[max_idx][i];
  }
  std::cout << "CHANGED BOUNDARY MAP TO AN EIGEN::VECTOR" << std::endl;
  Eigen::MatrixXd bdry_pos;
  igl::map_vertices_to_circle(V, bdry, bdry_pos);
  std::cout << "OK figured out what the boundary is" << std::endl;

  // Now apply the weighted Tutte map
  Eigen::MatrixXi E;
  igl::edges(F, E);
  std::cout << "OK figured out what the edges are too " << std::endl;

  Eigen::SparseMatrix<double> L;
  typedef Eigen::Triplet<double> T;
  std::vector<T> triples;
  triples.reserve(4 * E.rows());
  for (int i = 0; i < E.rows(); i++) {
	double weight = 1.0 / ((V.row(E(i,0)) - V.row(E(i,1))).norm());
        triples.push_back(T(E(i,0), E(i,1), weight));
        triples.push_back(T(E(i,1), E(i,0), weight));
	triples.push_back(T(E(i,0), E(i,0),-weight));
	triples.push_back(T(E(i,1), E(i,1),-weight));
  }  	  
  int n_v = F.maxCoeff() + 1;
  L.resize(n_v, n_v);
  L.setFromTriplets(triples.begin(), triples.end());
  std::cout << "figured out the Laplacian" << std::endl;

  // Now solve the weighted optimization problem subject to boundary constraints U_x^T L U_x + U_y^T L U_y
  igl::min_quad_with_fixed_data<double> data;

  // No linear constraints or linear term in optimization objective;
  Eigen::SparseMatrix<double> Z;
  Z.setZero();
  Eigen::VectorXd B;
  B.setZero();
  Eigen::VectorXd z(n_v);
  z.setZero();

  igl::min_quad_with_fixed_precompute(L, bdry, Z, false, data); // no linear equality constraints
  std::cout << "factorized data out laplacian" << std::endl;
  U.resize(n_v, 2);
  for (int i = 0; i <= 1; i++) {
	Eigen::VectorXd cur_col;
	igl::min_quad_with_fixed_solve(data, z, bdry_pos.col(i), B, cur_col); // no linear term in objective or linear inequality contraints
  	std::cout << "I SOLVED IT " << i << std::endl;
  	U.col(i) = cur_col;
   }
  std::cout << "figured out full paramterization" << std::endl;
}

