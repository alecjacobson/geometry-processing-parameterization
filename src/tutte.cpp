#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  /*
   We want to map the vertices onto a plane in order to parameterize the surface.
   Thinking of it as a spring system:
   
   summing over all of the edges (i,j)ϵE
   min ∑ ‖uᵢ - uⱼ‖²
   
   will give trivial solutions unless some points are constrained (in our case we will map the boundary edges to a unit circle), 
   and we also want to maintain the relative lengths of the edges when mapped:

   min ∑ wᵢⱼ‖uᵢ - uⱼ‖²

   our wᵢⱼ will come from the cotangent laplacian, and can rewrite the problem as:

   min ½ tr(UᵀLU)
   min ½ (U₁ᵀLU₁ + U₂ᵀLU₂)

   which we can solve indep with min_quad_with_fixed
   */

  U = Eigen::MatrixXd(V.rows(),2);
  Eigen::VectorXi boundaryEdges;
  igl::boundary_loop(F, boundaryEdges);

  Eigen::MatrixXd boundaryPositions;
  igl::map_vertices_to_circle(V, boundaryEdges, boundaryPositions);

  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // min_quad_with_fixed(A, B, known, Y, Aeq, Beq, pd, Z)
  // tries to minimize WRT Z:  ½ ZᵀAZ + BZ + C
  // with constraints: Aeq*Z = Beq
  Eigen::SparseMatrix<double> Aeq = Eigen::SparseMatrix<double>(0,0);
  Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(0,0);
  Eigen::VectorXd B = Eigen::VectorXd::Zero(L.rows());
  Eigen::VectorXd Uu, Uv;
  Eigen::VectorXd YOne = boundaryPositions.col(0);
  Eigen::VectorXd YTwo = boundaryPositions.col(1);
  // igl::min_quad_with_fixed(L, B, boundaryEdges, YOne, Aeq, Beq, false, Uu);
  // igl::min_quad_with_fixed(L, B, boundaryEdges, YTwo, Aeq, Beq, false, Uv);
  igl::min_quad_with_fixed_data<double> mqwf;
  igl::min_quad_with_fixed_precompute(L, boundaryEdges, Aeq, false, mqwf);
  igl::min_quad_with_fixed_solve(mqwf, B, YOne, Beq, Uu);
  igl::min_quad_with_fixed_solve(mqwf, B, YTwo, Beq, Uv);
  
  U.col(0) = Uu;
  U.col(1) = Uv;
}

