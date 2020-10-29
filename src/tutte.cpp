#include "tutte.h"

#include "igl/map_vertices_to_circle.h"
#include "igl/boundary_loop.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/cotmatrix.h"

void tutte(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &U) {


    // Let's solve for min(0.5*U'*L*U) such that
    // all boundary vertices should be mapped to the unit circle and
    // interior vertices mapped inside the disk without flips.


    // I am assuming wij=1. Computing wij based on edge-lengths is explained as an alternative(hack) but not requested in the assignment/implementation.
    // Therefore L is the uniform Laplacian and can be easily computed using the whitelisted functions.
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    // Overloaded function (boundary_loop).
    // let's pass a Eigen::PlainObjectBase" (eg. vectorXi) as second argument
    // so that it directly returns the longest loop in terms of vertices.
    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);


    // mapping bnd to a circle
    Eigen::MatrixXd UV;
    igl::map_vertices_to_circle(V, bnd, UV);

    // Optimization

    // igl::min_quad_with_fixed minimizes
    // trace( 0.5*Z'*A*Z + Z'*B + constant )
    // subject to
    // Z(known,:) = Y, and
    // Aeq*Z = Beq

    // This is super handy because
    // We want to solve trace( 0.5*U'*L*U + U'*0 + 0*constant )
    // = trace( 0.5*U'*L*U)
    // subject to
    // U(bnd,:)=UV -> points in the boundary attached to the circle.
    // so ...
    Eigen::VectorXd B = Eigen::VectorXd::Zero(V.rows());
    Eigen::SparseMatrix<double> constant; //this is init to zeros since it is an sp. matrix...
    igl::min_quad_with_fixed(L, B, bnd, UV, constant, B, false, U);
}

