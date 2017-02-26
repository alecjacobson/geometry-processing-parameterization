#include "tutte.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <igl/cotmatrix.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/slice_into.h>

#include <vector>
#include <ostream>

void tutte2(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & U );


void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    // Replace with your code
    U = V.leftCols(2); // projects onto xy plane with foldovers
    
    // min_U 1/2 tr( U^T L U )
    // Can use cotmatrix for L
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix( V, F, L );
    
    // can make w_ij = 1 / || v_i - v_j ||

    // constrain the boundary to the outside of the circle.
    //std::vector< std::vector<int> > boundaryVerts
    // Note - we are assuming there's only 1 boundary that matters, the longest.
    Eigen::VectorXi boundaryVerts;
    igl::boundary_loop( F, boundaryVerts );

    Eigen::MatrixXd boundaryUV;
    igl::map_vertices_to_circle( V, boundaryVerts, boundaryUV );


    //igl::min_quad_with_fixed( )
    igl::min_quad_with_fixed_data<double> mqwf;
    Eigen::VectorXd B = Eigen::VectorXd::Zero( V.rows(), 1 ); // linear term
    Eigen::VectorXd Beq; // empty constraints
    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute( (-L).eval(), boundaryVerts, Aeq, true, mqwf);

    Eigen::VectorXd U1 = U.col(0), U2 = U.col(1);
    igl::min_quad_with_fixed_solve( mqwf, B, boundaryUV.col(0), Beq, U1 );
    igl::min_quad_with_fixed_solve( mqwf, B, boundaryUV.col(1), Beq, U2 );

    U.col(0) = U1;
    U.col(1) = U2;

}

void tutte2(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & U)
{
/*    // cut from tutte....
    Eigen::VectorXi allVerts;
    igl::colon<int>( 0, V.rows()-1, allVerts );
    
    // make a list of verts to use:
    // in - inside or free variable verts
    // b - boundary or assigned verts (set on our outside circle)
    Eigen::VectorXi in, IA;
    igl::setdiff( allVerts, boundaryVerts, in, IA );
    
    // Construct and slice up Laplacian
    Eigen::SparseMatrix<double> L_in_in, L_in_b;
    igl::slice( L, in, in, L_in_in );
    igl::slice( L, in, boundaryVerts, L_in_b );
    
    // Solve PDE - for U1 (u)
    //Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver( -L_in_in );
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double > > solver( -L_in_in );
    Eigen::VectorXd u_in = solver.solve( L_in_b * boundaryUV.col(0) ); // L_in_b*bc
    // slice into solution
    Eigen::VectorXd U1 = U.col(0);
    igl::slice_into( u_in, in, U1 );

    Eigen::VectorXd v_in = solver.solve( L_in_b * boundaryUV.col(1) ); // L_in_b*bc
    // slice into solution
    Eigen::VectorXd U2 = U.col(1);
    igl::slice_into( v_in, in, U2 );

    U.col(0) = U1;
    U.col(1) = U2;

*/    
}
