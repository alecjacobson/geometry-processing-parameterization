#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <vector>
#include <set>
#include <iterator>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include <igl/map_vertices_to_circle.h>
#include "active_set.h"


void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  U = V.leftCols(2);

  U.setZero();

  using Index = Eigen::MatrixXi::Scalar;
  std::vector<Index> bnd_loop;
  igl::boundary_loop(F,bnd_loop);

  Eigen::VectorXi bl(bnd_loop.size());
  std::copy(bnd_loop.begin(),bnd_loop.end(),bl.data());
  Eigen::MatrixXd UV;

  igl::map_vertices_to_circle(V,bl,UV);


  ActiveSet coords(U.rows());

  coords.setVariableType<ActiveSet::VariableType::Dirichlet>(bnd_loop);
  coords.updateActive();







  Eigen::SparseMatrix<double> L,A,D;
  igl::cotmatrix(V,F,L);

  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,D);


  int active_size = coords.active_size();


  A.resize(active_size,active_size);

  U = coords.expand<ActiveSet::VariableType::Dirichlet>(UV);

  Eigen::MatrixXd pu = L * U;



  Eigen::MatrixXd rhs = coords.reduce<ActiveSet::VariableType::Free>(pu);


  A = coords.reduce<ActiveSet::VariableType::Free>(L);


  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
  solver.compute(A);
  for(int i = 0 ; i < 2; ++i) {
      UV.col(i) = solver.solve(rhs.col(i));
      std::cout << "RHS("<<i<<"): ";
      std::cout << rhs.col(i).transpose() << std::endl;
      std::cout << "SLN_rhs("<<i<<"): ";
      std::cout << (A * UV.col(i) - rhs.col(i)) << std::endl;

  }

  U += coords.expand<ActiveSet::VariableType::Free>(UV);




}


