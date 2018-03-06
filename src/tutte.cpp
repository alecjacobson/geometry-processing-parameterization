#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U) {

    // resize U to contain a row for each vertex in V and columns corresponding
    // to the (u,v) coordinates for that vertex
    U.resize(V.rows(), 2);
    
    // compute a list of boundary vertices on the longest boundary loop
    Eigen::VectorXi bound;
    igl::boundary_loop(F, bound);
    
    // map the boundary vertices to the unit circle
    Eigen::MatrixXd circ;
    igl::map_vertices_to_circle(V, bound, circ);
    
    // compute the cotangent Laplacian, to be used in the minimization problem
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);

    // min_quad_with_fixed_precompute populates data with the necessary information to solve
    // a minimization of a quadratic energy of the form: tr(0.5*Z'*A*Z + Z'*B + constant) 
    // (where ' denotes transpose) subject to Y constraints on Z and linear
    // equality constraints Aeq*Z = Beq. Since our minimization problem
    // is of the form 0.5*tr(U'LU), we get the following correspondences: L -> A; U -> Z;
    // B -> Zero. Likewise, Aeq and Beq are empty since we do not impose linear
    // equality constraints.
    Eigen::SparseMatrix<double> Aeq;
    Eigen::MatrixXd Beq;
    Eigen::VectorXd B = Eigen::VectorXd::Zero(L.rows());
    igl::min_quad_with_fixed_data<double> data;
    igl::min_quad_with_fixed_precompute(L, bound, Aeq, false, data);
    
    // solve the minimization problem
    igl::min_quad_with_fixed_solve(data, B, circ, Beq, U);
}

