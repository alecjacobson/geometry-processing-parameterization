#include "vector_area_matrix.h"



#include <igl/boundary_facets.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
    // Replace with your code
    //A.resize(F.maxCoeff()+1,F.maxCoeff()+1);

    // max number of vertices referenced by faces
    const int n = F.maxCoeff() + 1;
    
    Eigen::MatrixXi E; // boundary edges
    igl::boundary_facets( F, E );

    
    std::vector< Eigen::Triplet<double> > trips;

    const int nEdges = E.rows();
    for( int e=0; e < nEdges; ++e )
    {
        int first  = E( e, 0 );
        int second = E( e, 1 );

        //  0   0   0  .5
        //  0   0 -.5   0
        //  0 -.5   0   0
        // .5   0   0   0
        trips.push_back( Eigen::Triplet<double>( first,      second + n,  .5 ) );
        trips.push_back( Eigen::Triplet<double>( first + n,  second,     -.5 ) );
        trips.push_back( Eigen::Triplet<double>( second,     first + n,  -.5 ) );
        trips.push_back( Eigen::Triplet<double>( second + n, first,       .5 ) );
    }

    Eigen::SparseMatrix< double > mtx( 2*n, 2*n );
    mtx.setFromTriplets( trips.begin(), trips.end() );

    A = mtx * 0.5;
}

//
// 
// Try and generate our A in a symmetric fashion...
// Our determinant for a pair of vertices on a boundary edge has a
// quartenary quadratic form.  If our 2 vertices are:
// U_1u, U_1v and  U_2u, U_2v
// we can stack these as {U_1u, U_1v, U_2u, U_2v}
// To write out our quartenary form, we can say that
// our points are { x, y, z, w }
// We'll use the conventional a,b,c... etc to be the coeffients.
// So our general quartenary quadratic would be:
//
// ax^2 + by^2 +cz^2 + dw^2 + exy + fxz + gxw + hyz +iyw + jzw
//
// When we write it all out, we get a matrix A for a pair of points
// that looks like
//
// 2a  e  f  g
// e  2b  h  i
// f   h 2c  j
// g   i  j 2d 
//
// Our determinant is: U_1u U_2v - U_2u U_1v
// which correspond to:
// xw - yz
// or -> gxw - hyz
//
// so in our matrix we need g=1, h=-1
// so A in this case becomes (with 1's across the diag instead of 2):
//
//  0   0   0  .5
//  0   0 -.5   0
//  0 -.5   0   0
// .5   0   0   0
//
// which we need to apply to each set of boundary edges.
// So we need to put this matrix into A for each set of  edges.
//
// In a way, this looks like two 90 degree vector rotations
//   0  .5
// -.5  0
//
//  and
//                             
//  0 -.5
// .5  0
//
// So we do the first rotation, then the second for the boundary points.
//   
