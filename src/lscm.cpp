#include "lscm.h"
#include "vector_area_matrix.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/repdiag.h"
#include "igl/eigs.h"
#include <Eigen/SVD>

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    int numV = F.maxCoeff() + 1;
    Eigen::SparseMatrix<double> A, Q;
    Eigen::MatrixXd EigVecs, tempU;
    Eigen::VectorXd EigVals;
    
    
    vector_area_matrix(F,A);
    
    Eigen::SparseMatrix<double> L, L2,M,B;
    L.resize(numV,numV);
    M.resize(numV,numV);
    
    L2.resize(2*numV,2*numV);
    Q.resize(2*numV,2*numV);
    B.resize(2*numV,2*numV);
    
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
    igl::repdiag(L,2,L2);
    igl::repdiag(M,2,B);

    Q = L2 - A;
    
    igl::eigs(Q,B,3,igl::EIGS_TYPE_SM, EigVecs, EigVals);
    
    cout << EigVals << "\n";
    cout << EigVecs.col(2).transpose() * B * EigVecs.col(2) << "\n";
    tempU.resize(numV,2);
    tempU.col(0) = EigVecs.col(2).head(numV);
    tempU.col(1) = EigVecs.col(2).tail(numV);
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svdOfU(tempU, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    U = svdOfU.matrixU();
    
    U.col(0) = U.col(0) * svdOfU.singularValues().row(0);
    U.col(1) = U.col(1) * svdOfU.singularValues().row(1);
    
}
