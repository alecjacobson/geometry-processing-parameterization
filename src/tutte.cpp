#include "tutte.h"
#include "igl/boundary_loop.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/map_vertices_to_circle.h"

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    
    Eigen::VectorXi boundary;
    Eigen::MatrixXd Y,Beq, B;
    Eigen::SparseMatrix<double> A, Aeq;
    
  //Need to parametrize the loop
    igl::boundary_loop(F, boundary);
    
    
    int numB = boundary.rows();
    
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(F.rows()*12 + numB*4);
    
    Eigen::MatrixXd vals(3,2);
    
    //Makes it easier to reference other vertices
    vals(0,0) = 1;
    vals(0,1) = 2;
    vals(1,0) = 0;
    vals(1,1) = 2;
    vals(2,0) = 0;
    vals(2,1) = 1;
    
    int numV = F.maxCoeff() + 1;
    
    A.resize(numV,numV);
    
    for (int i = 0; i < F.rows(); i ++) {
        for (int j = 0; j < 3; j ++) {
            //coeffs for u coordinate
            tripletList.push_back(T(F(i,vals(j,0)), F(i,vals(j,1)), -1.0));
            tripletList.push_back(T(F(i,vals(j,1)), F(i,vals(j,0)), -1));
            
            tripletList.push_back(T(F(i,vals(j,0)), F(i,vals(j,0)), 1));
            tripletList.push_back(T(F(i,vals(j,1)), F(i,vals(j,1)), 1));
            
            
        }
        
        
    }
    
    int Val;
    for (int i = 0; i < numB; i ++) {
        if (i+1 == numB) {
            Val = 0;
        }
        else {
            Val = i + 1;
        }
            tripletList.push_back(T(boundary(i), boundary(Val), -1));
            tripletList.push_back(T(boundary(Val), boundary(i), -1));
            tripletList.push_back(T(boundary(i), boundary(i), 1));
            tripletList.push_back(T(boundary(Val), boundary(Val), 1));
        
    }

   
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    B = Eigen::MatrixXd::Zero(numV,2);
    
    Y.resize(numB,2);
    igl::map_vertices_to_circle(V,boundary, Y);
    
    Aeq.resize(1,numV);
    Beq = Eigen::MatrixXd::Zero(1,2);
    U.resize(numV,2);
    
    igl::min_quad_with_fixed(A,B, boundary, Y, Aeq,Beq, 1, U);
}

