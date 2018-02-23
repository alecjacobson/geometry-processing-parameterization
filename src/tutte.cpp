#include "tutte.h"
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <Eigen/src/IterativeLinearSolvers/ConjugateGradient.h>
#include <Eigen\src\Jacobi\Jacobi.h>

#include <iostream>

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
	int numRows = matrix.rows() - 1;
	int numCols = matrix.cols();

	if (rowToRemove < numRows)
		matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = 
			matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

	matrix.conservativeResize(numRows, numCols);
}


void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{

  Eigen::VectorXi bnd;
  igl::boundary_loop(F, bnd);

  Eigen::MatrixXi F_copy = F;
  Eigen::MatrixXd V_copy = V;

  Eigen::MatrixXd int_V = V;
  Eigen::MatrixXd split_V(V);

  Eigen::VectorXi map(V.rows());

  for (int i = 0; i < bnd.rows(); i++)
  {
	  split_V.row(i) = V.row(bnd(i));
	  map(i) = bnd(i);
	  removeRow(int_V, bnd(i));
	  V_copy.row(bnd(i)).setZero();
  }

  std::cout << "# of vertices on the boundary " << bnd.rows() << std::endl;
  std::cout << "# of vertices on the interior " << int_V.rows() << std::endl;
  std::cout << "# of vertices total " << V.rows() << std::endl;

  int index = 0;
  for (int i = 0; i < V.rows(); i++)
  {
	  if (!(V_copy.row(i).isZero()))
	  {
		  map(index + bnd.rows()) = i;
		  index++;
	  }
  }

  std::cout << "checking map size " << map.rows() << std::endl;

  for (int i = 0; i < F.rows(); i++)
  {
	  auto f = F_copy.row(i);
	  f(0) = map(f(0));
	  f(1) = map(f(1));
	  f(2) = map(f(2));
  }

  std::cout << "# of faces in the new matrix " << F_copy.rows() << std::endl;

  int bnd_size = bnd.rows();
  int int_verts_size = int_V.rows();
  split_V.block(bnd_size, 0, int_verts_size, 3) = int_V.block(0, 0, int_verts_size, 3);

  std::cout << "# of vertices in the split matrix " << split_V.rows() << std::endl;

  Eigen::SparseMatrix<double> laplacian(split_V.rows(),split_V.rows());
  igl::cotmatrix(split_V, F_copy, laplacian);

  std::cout << "size of cotan laplacian " << laplacian.rows() << ", " << laplacian.cols() << std::endl;

  Eigen::MatrixXd bnd_mapping;
  igl::map_vertices_to_circle(V, bnd, bnd_mapping);

  std::cout << "size of boundary mapping " << bnd_mapping.rows() << std::endl;

  Eigen::MatrixXd mapping(V.rows(), 2);
  U.resize(V.rows(), 2);
  mapping.block(0, 0, bnd_mapping.rows(), 2) = bnd_mapping;

  std::cout << "size of U " << mapping.rows() << std::endl;

  Eigen::MatrixXd int_laplacian = laplacian.block(bnd_size, bnd_size, int_verts_size, int_verts_size);
  std::cout << "size of  interior laplacian " << int_laplacian.rows() << std::endl;

  Eigen::MatrixXd A = laplacian.block(bnd_size, 0, int_verts_size, int_verts_size) * bnd_mapping;
  std::cout << bnd_size + int_verts_size << std::endl;
  Eigen::MatrixXd int_mapping;// = -int_laplacian.inverse() * laplacian.block(bnd_size, 0, int_verts_size, int_verts_size) * bnd_mapping;
  std::cout << "Ax=b "  << std::endl;

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
  //cg.compute(A);
  int_mapping = svd.solve(int_laplacian);
  std::cout << "did it work? " << std::endl;


  U.block(bnd_mapping.rows(), bnd_mapping.rows(), int_mapping.rows(), 2) = int_mapping;
  
  for (int i = 0; i < U.rows(); i++)
  {
	  U.row(i) = mapping.row(map(i));
  }
  
}

