#include "vector_area_matrix.h"
#include <Eigen/Eigen>
#include <igl/boundary_loop.h>
#include <igl/vector_area_matrix.h>
typedef Eigen::Triplet<double> T;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  int V_size = F.maxCoeff()+1;
  int n = V_size;
  A.resize(V_size*2,V_size*2);

  // Get boundary loop vertices:
  std::vector<std::vector<int>> bound_loop;
  igl::boundary_loop(F, bound_loop);

  // Calculate area matrix:
  int loop_num = bound_loop.size();
  std::vector<T> triplet_list;

  for (int i = 0; i < loop_num; i ++) {
    std::vector<int> cur_loop = bound_loop[i];
    int n = cur_loop.size();
    for (int j = 0; j < cur_loop.size(); j ++) {
      int ui = cur_loop[j];
      int uj = cur_loop[(j + 1) % n];

      // det(ui, uj):
      // uiu * ujv - uiv * uju
      triplet_list.push_back(T(ui, uj + V_size, 0.5));
      triplet_list.push_back(T(ui + V_size, uj, -0.5));
    }
  }

  // Symmetric:
  A.setFromTriplets(triplet_list.begin(), triplet_list.end());
  Eigen::SparseMatrix<double> AT =  A.transpose();
  A = A + AT; // times 0.5 here seems incorrect
}

