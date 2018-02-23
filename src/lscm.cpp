#include "lscm.h"

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    // Replace with your code
    int n = F.maxCoeff()+1;
    Eigen::SparseMatrix<double> Area;
    vector_area_matrix(F,Area);
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);
    Eigen::SparseMatrix<double> Q;
    igl::repdiag(L,2,Q);
    Q = Q-2*Area;
    
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    Eigen::SparseMatrix<double> B;
    igl::repdiag(M,2,B);
    
    //Adapted from the libigl tutorial
    Eigen::MatrixXd Vectors;
    Eigen::VectorXd Values;
    igl::eigs(Q,B,3,igl::EIGS_TYPE_SM,Vectors,Values);
    
    
    U.resize(n,2);
    U.col(0) = Vectors.col(2).head(n);
    U.col(1) = Vectors.col(2).tail(n);
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(U, Eigen::ComputeThinU | Eigen::ComputeThinV);
    U = svd.matrixU() * svd.singularValues().asDiagonal();
}
