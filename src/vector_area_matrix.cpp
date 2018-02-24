#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  
  int n = F.maxCoeff() + 1;
  A.resize(2*n, 2*n);

  std::vector<std::vector<int>> bnds;
  igl::boundary_loop(F, bnds);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(n);

  for (int i = 0; i < bnds.size(); i++)
  {
	  std::vector<int> bnd = bnds[i];
	  for (int j = 0; j < bnd.size(); j++)
	  {
		  //std::cout << j << std::endl;

		  int index_1 = bnd[j];
		  int index_2 = 0;
		  if (j == bnd.size()-1)
			  index_2 = bnd[0];
		  else
			  index_2 = bnd[j + 1];
		  
		  //if (index_1+n > n*2 || index_2+n > n*2)
			//std::cout << index_1+n << "HEYHEYHYEHEYHEYHEY " << index_2+n << std::endl;
		  tripletList.push_back(T(index_1, index_2+n, 1));
		  tripletList.push_back(T(index_2, index_1+n, -1));
	  }
  }

  A.setFromTriplets(tripletList.begin(), tripletList.end());

  Eigen::SparseMatrix<double> B;
  B = Eigen::SparseMatrix<double>(A.transpose()) + A;
  
  A = 0.5 * B;
  
}

