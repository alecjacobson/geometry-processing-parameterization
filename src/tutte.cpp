#include "tutte.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/edges.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = V.leftCols(2);
  Eigen::VectorXi bnd;
  igl::boundary_loop(F,bnd);

  Eigen::MatrixXi E;
  igl::edges(F,E);
  typedef Eigen::Triplet<double> T;
  std::vector<T> tlist;
  tlist.reserve(E.rows()*2);
  for (int i=0; i<E.rows();i++) {
    tlist.push_back(T(E(i,0),E(i,1),1.0));
    tlist.push_back(T(E(i,1),E(i,0),1.0));
  }
  Eigen::SparseMatrix<double> L;
  int L_size = E.maxCoeff() + 1;
  L.resize(L_size,L_size);
  L.setFromTriplets(tlist.begin(),tlist.end());

  for(int i=0; i<L_size; i++) {
    L.insert(i,i) = -L.row(i).sum(); 
  }

  Eigen::MatrixXd UV;
  igl::map_vertices_to_circle(V,bnd,UV);

  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(V.rows(),2);
  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(L,bnd,Eigen::SparseMatrix<double>(),false,data);
  igl::min_quad_with_fixed_solve(data,B,UV,Eigen::MatrixXd(),U);
}