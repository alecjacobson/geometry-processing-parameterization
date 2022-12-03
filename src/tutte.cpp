#include "tutte.h"
#include "igl/boundary_loop.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/edges.h"
#include "igl/min_quad_with_fixed.h"

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Determine the longest boundary loop and map its vertices onto the unit circle
  Eigen::VectorXi bdry;
  igl::boundary_loop(F, bdry);
  Eigen::MatrixXd bdry_pos;
  igl::map_vertices_to_circle(V, bdry, bdry_pos);

  // Calculate the weighted Laplacian matrix
  Eigen::MatrixXi E;
  igl::edges(F, E);
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

  // Now solve the weighted optimization problem  U_x^T L U_x + U_y^T L U_y subject to boundary constraints calculated earlier
  igl::min_quad_with_fixed_data<double> data;

  // No linear constraints or linear term in optimization objective;
  Eigen::SparseMatrix<double> Z;
  Z.setZero();
  Eigen::VectorXd B;
  B.setZero();
  Eigen::VectorXd z(n_v);
  z.setZero();

  igl::min_quad_with_fixed_precompute(L, bdry, Z, false, data); // no linear equality constraints
  U.resize(n_v, 2);
  for (int i = 0; i <= 1; i++) {
	Eigen::VectorXd cur_col;
	igl::min_quad_with_fixed_solve(data, z, bdry_pos.col(i), B, cur_col); // no linear term in objective or linear inequality contraints
  	U.col(i) = cur_col;
   }
}

