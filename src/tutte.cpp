#include "tutte.h"
#include "math.h"
#include <vector>
#include <igl/boundary_loop.h>
#include <igl/edges.h>
#include <fstream>
#include <igl/min_quad_with_fixed.h>
#include <igl/repdiag.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  std::vector<int> L;
  igl::boundary_loop(F,L);
  Eigen::MatrixXi E;
  igl::edges(F,E);
  double tp;
  typedef Eigen::Triplet<double> T;
  std::vector<T> list;
  Eigen::VectorXd diag=Eigen::VectorXd::Zero(V.rows());
  Eigen::VectorXi b(L.size());
  double pi=3.14159265359;
  Eigen::MatrixXd bc(L.size(),2);
  for (int i=0; i<L.size(); i++){
    bc.row(i)=Eigen::RowVector2d(cos(pi*2/L.size()*i),sin(pi*2/L.size()*i));
    b(i)=L[i];
  }
  for (int i=0; i<E.rows(); i++){
  	tp=1.0/(V.row(E(i,0))-V.row(E(i,1))).norm();
  	list.push_back(T(E(i,0),E(i,1),tp));
  	list.push_back(T(E(i,1),E(i,0),tp));
  	diag(E(i,0))-=tp;
  	diag(E(i,1))-=tp;
  }
  for (int i=0; i<V.rows(); i++)
    list.push_back(T(i,i,diag(i)));
  Eigen::SparseMatrix<double> A;
  A.resize(V.rows(),V.rows());
  A.setFromTriplets(list.begin(), list.end()); 
  U.resize(V.rows(),2);
  const Eigen::VectorXd B_flat = Eigen::VectorXd::Zero(V.rows());
  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(A,b,Eigen::SparseMatrix<double>(),false,data);
  igl::min_quad_with_fixed_solve(data,B_flat,bc,Eigen::VectorXd(),U);
}

