#include "tutte.h"
#include "igl/edge_lengths.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/boundary_loop.h"

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);
  Eigen::SparseMatrix<double> L(V.rows(), V.rows());
  for (int i = 0; i < F.rows(); i++) {
  	//[1,2],[2,0],[0,1]
  	L.coeffRef(F(i, 0), F(i, 1)) += 1 / l(i, 2);
  	L.coeffRef(F(i, 1), F(i, 2)) += 1 / l(i, 0);
  	L.coeffRef(F(i, 2), F(i, 0)) += 1 / l(i, 1);

    L.coeffRef(F(i, 1), F(i, 0)) = L.coeffRef(F(i, 0), F(i, 1));
    L.coeffRef(F(i, 2), F(i, 1)) = L.coeffRef(F(i, 1), F(i, 2));
    L.coeffRef(F(i, 0), F(i, 2)) = L.coeffRef(F(i, 2), F(i, 0));

    L.coeffRef(F(i, 0), F(i, 0)) -= 1 / l(i, 2);
    L.coeffRef(F(i, 1), F(i, 1)) -= 1 / l(i, 2);
    L.coeffRef(F(i, 1), F(i, 1)) -= 1 / l(i, 0);
    L.coeffRef(F(i, 2), F(i, 2)) -= 1 / l(i, 0);
    L.coeffRef(F(i, 2), F(i, 2)) -= 1 / l(i, 1);
    L.coeffRef(F(i, 0), F(i, 0)) -= 1 / l(i, 1);
  }
  Eigen::VectorXi knowns;
  igl::boundary_loop(F, knowns);
  Eigen::MatrixXd bc;
  igl::map_vertices_to_circle(V, knowns, bc);
  igl::min_quad_with_fixed_data<double> mqwf;
  Eigen::SparseMatrix<double> Aeq(0, 0);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(V.rows(), 1);
  Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(0, 0);
  igl::min_quad_with_fixed_precompute(L, knowns, Aeq, false, mqwf);
  Eigen::VectorXd U0 = Eigen::VectorXd::Zero(V.rows());
  Eigen::VectorXd U1 = Eigen::VectorXd::Zero(V.rows());
  igl::min_quad_with_fixed_solve(mqwf,B,bc.col(0),Beq,U0);
  igl::min_quad_with_fixed_solve(mqwf,B,bc.col(1),Beq,U1);
  U = Eigen::MatrixXd::Zero(V.rows(), 2);
  U.col(0) = U0;
  U.col(1) = U1;
}

