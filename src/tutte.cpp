#include "tutte.h"

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  int n = F.maxCoeff() + 1;

  Eigen::MatrixXi edges;
  igl::edges(F, edges);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(edges.rows()*4);

  // construst weighted graph laplacian
  Eigen::SparseMatrix<double> L;
  for (int i = 0; i < edges.rows(); i++){
    // loop over edges
    int pi = edges(i,0);
    int pj = edges(i,1);
    double wij = 1.0 / (V.row(pi) - V.row(pj)).norm();

    tripletList.push_back(T(pi, pj, wij));
    tripletList.push_back(T(pj, pi, wij));
    tripletList.push_back(T(pi, pi, -wij));
    tripletList.push_back(T(pj, pj, -wij));
  }
  L.resize(n, n);
  L.setFromTriplets(tripletList.begin(), tripletList.end());
  
  Eigen::VectorXi bdry;
  // map_vertices_to_circle takes Eigen::VectorXi for bdry
  igl::boundary_loop(F, bdry);

  Eigen::MatrixXd C; //constraints
  igl::map_vertices_to_circle(V, bdry, C);

  Eigen::VectorXd B = Eigen::VectorXd::Zero(n);
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed(L, B, bdry, C, Aeq, B, false, U);
}
