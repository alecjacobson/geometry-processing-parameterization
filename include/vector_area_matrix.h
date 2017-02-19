#ifndef VECTOR_AREA_MATRIX_H
#define VECTOR_AREA_MATRIX_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
// Constructs the symmetric area matrix A, s.t.  [V.col(0)' V.col(1)'] * A *
// [V.col(0); V.col(1)] is the **vector area** of the mesh (V,F).
//
// Inputs:
//   F  #F by 3 list of mesh faces (must be triangles)
// Outputs:
//   A  #Vx2 by #Vx2 sparse area matrix
//
void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A);

#endif

