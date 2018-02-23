#include "tutte.h"

void graph_laplacian(
                     const Eigen::MatrixXi & E,
                     Eigen::SparseMatrix<double> & L)
{
    int n = E.maxCoeff() + 1;
    L.resize(n,n);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(E.rows()*2);
    
    for(int e = 0; e < E.rows(); e++){
        tripletList.push_back(T(E(e,0),E(e,1),1.0));
        tripletList.push_back(T(E(e,1),E(e,0),1.0));
        tripletList.push_back(T(E(e,1),E(e,1),-1.0));
        tripletList.push_back(T(E(e,0),E(e,0),-1.0));
    }
    
    L.setFromTriplets(tripletList.begin(), tripletList.end());
}

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    // Replace with your code
    Eigen::VectorXi Loop;
    igl::boundary_loop(F,Loop);
    
    Eigen::MatrixXd UV;
    igl::map_vertices_to_circle(V,Loop,UV);
    
    //Get cot matrix
    Eigen::SparseMatrix<double> L;
    Eigen::MatrixXi E;
    igl::edges(F,E);
    graph_laplacian(E,L);
//    igl::cotmatrix(V,F,L);
    
    
    U = Eigen::MatrixXd::Zero(V.rows(),2);
    //Adapted from the libigl tutorial:
    Eigen::VectorXd B = Eigen::VectorXd::Zero(V.rows());
    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_data<double> mqwf;
    igl::min_quad_with_fixed_precompute(
                                    L,Loop,Aeq,false,mqwf);
    Eigen::VectorXd Beq;
    igl::min_quad_with_fixed_solve<double>(mqwf,B,UV,Beq,U);
}

