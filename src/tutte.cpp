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
    
    //Finds elements on the boundary
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
    
    //We define A in two steps
    //We loop through every face and add a triplet that will add the correct terms to the correct places.
    //This will double count for edges not on the boundary, so we have another loop to double count edges on the boundary
    //This way we avoid have to do any collision detection or something to ensure no over or under counting occurs
    
    //Looping over all the edges and decomposing (a-b)^2 = a^2 - ab - ba + b^2 to fill in cells sparselu
    for (int i = 0; i < F.rows(); i ++) {
        for (int j = 0; j < 3; j ++) {
            
            //off diagonal entries
            tripletList.push_back(T(F(i,vals(j,0)), F(i,vals(j,1)), -1.0));
            tripletList.push_back(T(F(i,vals(j,1)), F(i,vals(j,0)), -1));
            
            //diagonal entries
            tripletList.push_back(T(F(i,vals(j,0)), F(i,vals(j,0)), 1));
            tripletList.push_back(T(F(i,vals(j,1)), F(i,vals(j,1)), 1));
            
            
        }
        
        
    }
    
    //Reading the boundary terms because they are only encountered once
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

   //Generate triplet
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    
    //Set linear term to be 0
    B = Eigen::MatrixXd::Zero(numV,2);
    
    //Map boundary terms to circle
    Y.resize(numB,2);
    igl::map_vertices_to_circle(V,boundary, Y);
    
    //Set trivial matrices
    Aeq.resize(1,numV);
    Beq = Eigen::MatrixXd::Zero(1,2);
    U.resize(numV,2);
    
    //Run min quad 
    igl::min_quad_with_fixed(A,B, boundary, Y, Aeq,Beq, 1, U);
}

