#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <vector>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <igl/map_vertices_to_circle.h>
#include "active_set.h"


/*
#include <chrono>
struct timer {

    std::chrono::time_point<std::chrono::system_clock> 
        start = std::chrono::system_clock::now();
    ~timer() {
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> diff = end-start;
      std::cout << "duration: " << diff.count() << "s\n";
  }
};
*/

//#define USE_ALEC_SOLVER
#ifdef USE_ALEC_SOLVER
#include <igl/min_quad_with_fixed.h>
#endif


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




  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);

  Eigen::SparseMatrix<double> A;


#ifdef USE_ALEC_SOLVER
  igl::cotmatrix(V,F,L);

  A = -L ;


  igl::min_quad_with_fixed_data<double> data;
  min_quad_with_fixed_precompute(A,bl,Eigen::SparseMatrix<double>(),true,data);

  Eigen::VectorXd B = Eigen::VectorXd::Zero(U.rows(),1);
  for(int i = 0 ; i < 2; ++i) {
    const Eigen::VectorXd c = UV.col(i);
    Eigen::VectorXd u;
    min_quad_with_fixed_solve(data,B,c,Eigen::VectorXd(),u);
    U.col(i) = u;
  }



#endif
#ifndef USE_ALEC_SOLVER

   
  ActiveSet coords(U.rows());

  coords.setVariableType<ActiveSet::VariableType::Dirichlet>(bnd_loop);
  coords.updateActive();





  int active_size = coords.active_size();


  A.resize(active_size,active_size);

  U = coords.expand<ActiveSet::VariableType::Dirichlet>(UV);




  Eigen::MatrixXd rhs;
  
  rhs = coords.reduce<ActiveSet::VariableType::Free>((-L*U).eval());


  UV.resize(active_size,2);
  A = coords.reduce<ActiveSet::VariableType::Free>(L);


  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

  solver.compute(A);
  for(int i = 0 ; i < 2; ++i) {
      UV.col(i) = solver.solve(rhs.col(i));

  }

  U += coords.expand<ActiveSet::VariableType::Free>(UV);




#endif
}


