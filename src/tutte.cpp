#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/map_vertices_to_circle.h>

typedef Eigen::Triplet<double> T;

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  //U = V.leftCols(2);
  int v = F.maxCoeff()+1;

  //B in vector form for mapping vertices to circle
  Eigen::VectorXi B;
  igl::boundary_loop(F,B);
  std::vector<T> tlist;

  Eigen::SparseMatrix<double> L, Aeq;
  L.resize(v,v);

  //uniform laplacian
  for (int i = 0; i < F.rows(); i++){
    int v1 = F(i,0);
    int v2 = F(i,1);
    int v3 = F(i,2);

    //adjacent -1
    tlist.push_back(T(v1,v2,-1));
    tlist.push_back(T(v2,v1,-1));
    tlist.push_back(T(v1,v3,-1));
    tlist.push_back(T(v3,v1,-1));
    tlist.push_back(T(v2,v3,-1));
    tlist.push_back(T(v3,v2,-1));

    //diagonal degree 2
    tlist.push_back(T(v1,v1, 2));
    tlist.push_back(T(v2,v2, 2));
    tlist.push_back(T(v3,v3, 2));
  }

  L.setFromTriplets(tlist.begin(),tlist.end());

  Eigen::MatrixXd C;
  igl::map_vertices_to_circle(V,B,C);

  //solve U
  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(L, B, Aeq, false, data);

  Eigen::MatrixXd Beq, tempB;
  tempB.setZero(v,2);
  igl::min_quad_with_fixed_solve(data, tempB, C, Beq, U);
}

